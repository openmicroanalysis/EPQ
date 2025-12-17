package gov.nist.microanalysis.EPQTools;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.stream.Collectors;

import gov.nist.microanalysis.EPQLibrary.AtomicShell;
import gov.nist.microanalysis.EPQLibrary.EPQException;
import gov.nist.microanalysis.EPQLibrary.Element;
import gov.nist.microanalysis.EPQLibrary.ToSI;
import gov.nist.microanalysis.EPQLibrary.TransitionEnergy;
import gov.nist.microanalysis.EPQLibrary.XRayTransition;
import gov.nist.microanalysis.EPQLibrary.XRayTransitionSet;

/**
 * <p>
 * A small class for representing KLM lines.
 * </p>
 * <p>
 * Copyright: Pursuant to title 17 Section 105 of the United States Code this
 * software is not subject to copyright protection and is in the public domain
 * </p>
 * <p>
 * Company: National Institute of Standards and Technology
 * </p>
 *
 * @author Daniel "Ooblioob" Davis, Nicholas W. M. Ritchie
 * @version 1.0
 */

abstract public class KLMLine implements Comparable<KLMLine> {

   enum LabelType {
      ELEMENT("Element"), //
      ELEMENT_ABBREV("Element abbreviation"), //
      LARGE_ELEMENT("Large abbreviation"), //
      SIEGBAHN("Siegbahn"), //
      IUPAC("IUPAC"), //
      FAMILY("Family"), //
      NONE("No labels");

      private final String mName;

      private LabelType(String id) {
         mName = id;
      }

      @Override
      public String toString() {
         return mName;
      }
   }

   enum KLMLineType {
      InvalidType, Satellite, SumPeak, EscapePeak, KEdge, LEdge, MEdge, NEdge, KTransition, LTransition, MTransition, NTransition;

      public boolean isTransition() {
         return (this == KTransition) || (this == LTransition) || (this == MTransition) || (this == NTransition);
      }
   }

   /**
    * Returns a shell with which to associate this KLMLine
    */
   abstract public AtomicShell getShell();

   abstract public String toLabel(LabelType lt);

   protected double mEnergy;
   protected double mAmplitude;

   public static class Transition extends KLMLine implements Comparable<KLMLine> {
      private final XRayTransition mTransition;

      public Transition(XRayTransition xrt) throws EPQException {
         super(TransitionEnergy.Default.compute(xrt), xrt.getWeight(XRayTransition.NormalizeKLM));
         mTransition = xrt;
      }

      @Override
      public String toLabel(LabelType lt) {
         switch (lt) {
            case ELEMENT :
               return mTransition.getElement().toString();
            case LARGE_ELEMENT :
            case ELEMENT_ABBREV :
               return mTransition.getElement().toAbbrev();
            case SIEGBAHN :
               return mTransition.getSiegbahnName();
            case IUPAC :
               return mTransition.getIUPACName();
            case FAMILY :
               return mTransition.getElement().toAbbrev() + " " + XRayTransitionSet.getClassNameForTransition(mTransition.getTransitionIndex());
            case NONE :
            default :
               return "";
         }
      }

      @Override
      public String toString() {
         return mTransition.toString();
      }

      @Override
      public AtomicShell getShell() {
         return mTransition.getDestination();
      }

      public XRayTransition getTransition() {
         return mTransition;
      }

      @Override
      public boolean contains(Element elm) {
         return mTransition.getElement().equals(elm);
      }

      @Override
      public String toSiegbahn() {
         return mTransition.getSiegbahnName();
      }

      @Override
      public KLMLineType getType() {
         KLMLine.KLMLineType res = KLMLine.KLMLineType.InvalidType;
         switch (getShell().getFamily()) {
            case AtomicShell.KFamily :
               res = KLMLine.KLMLineType.KTransition;
               break;
            case AtomicShell.LFamily :
               res = KLMLine.KLMLineType.LTransition;
               break;
            case AtomicShell.MFamily :
               res = KLMLine.KLMLineType.MTransition;
               break;
            case AtomicShell.NFamily :
               res = KLMLine.KLMLineType.NTransition;
               break;
            default :
               res = KLMLine.KLMLineType.InvalidType;
               break;
         }
         return res;
      }

      /**
       * compareTo - Allows to KLM lines to be ordered.
       *
       * @param kl
       *           KLMLine
       * @return int
       */
      @Override
      public int compareTo(KLMLine kl) {
         // Ordering: Transition, SumPeak, EscapePeak, Edge
         if (kl instanceof Transition tr) {
            return this.mTransition.compareTo(tr.mTransition);
         } else {
            return 1; // Transitions are always come before other types
            // return (kl instanceof Transition) || (kl instanceof SumPeak) ||
            // (kl instanceof EscapePeak) || (kl instanceof Edge) ? -1 : 1;
         }
      }

      @Override
      public boolean equals(Object obj) {
         if (obj == null) {
            return false;
         }
         if (this == obj) {
            return true;
         }
         if (obj instanceof Transition tr) {
            return this.mTransition.equals(tr.mTransition);
         }
         return false;
      }
   }

   public static Set<Transition> suggestKLM(Element elm, double eMax) {
      final int[] DEFAULT_KLM_LINES = {XRayTransition.KA1, XRayTransition.KB1, XRayTransition.LA1, XRayTransition.LB1, XRayTransition.MA1,
            XRayTransition.MA2, XRayTransition.MB, XRayTransition.MG};
      final Set<Transition> lines = new TreeSet<>();
      for (final int tr : DEFAULT_KLM_LINES) {
         try {
            if (XRayTransition.exists(elm, tr) && (XRayTransition.getEnergy(elm, tr) < eMax)) {
               lines.add(new KLMLine.Transition(new XRayTransition(elm, tr)));
            }
         } catch (EPQException e) {
            // Ignore
         }
      }
      return lines;
   }

   public static class Edge extends KLMLine implements Comparable<KLMLine> {

      private final AtomicShell mShell;

      static public AtomicShell mostOccupiedShellInFamily(AtomicShell sh) {
         AtomicShell res = sh;
         try {
            final int fam = sh.getFamily();
            for (int shell = AtomicShell.getFirstInFamily(fam); shell <= AtomicShell.getLastInFamily(fam); shell++) {
               final AtomicShell tmp = new AtomicShell(sh.getElement(), shell);
               if (tmp.getGroundStateOccupancy() > res.getGroundStateOccupancy()) {
                  res = tmp;
               }
            }
         } catch (final Exception e) {
            e.printStackTrace();
         }
         return res;
      }

      static public double fractionalOccupancy(AtomicShell shell) {
         final double den = 1.0 / mostOccupiedShellInFamily(shell).getGroundStateOccupancy();
         return Double.isNaN(den) ? 1.0 : den * shell.getGroundStateOccupancy();
      }

      public Edge(AtomicShell shell) {
         super(shell.getEdgeEnergy(), 0.2 * fractionalOccupancy(shell));
         mShell = shell;
      }

      @Override
      public String toLabel(LabelType lt) {
         switch (lt) {
            case ELEMENT :
               return mShell.getElement().toString() + " edge";
            case LARGE_ELEMENT :
            case ELEMENT_ABBREV :
               return mShell.getElement().toAbbrev() + " edge";
            case SIEGBAHN :
               return mShell.getSiegbahnName() + " edge";
            case IUPAC :
               return mShell.getIUPACName() + " edge";
            case FAMILY :
               return mShell.getIUPACName() + " edge";
            case NONE :
            default :
               return "";
         }
      }

      @Override
      public AtomicShell getShell() {
         return mShell;
      }

      @Override
      public String toString() {
         return mShell.toString();
      }

      @Override
      public boolean contains(Element elm) {
         return mShell.getElement().equals(elm);
      }

      @Override
      public String toSiegbahn() {
         return mShell.toString();
      }

      @Override
      public KLMLineType getType() {
         KLMLine.KLMLineType res = KLMLine.KLMLineType.InvalidType;
         switch (getShell().getFamily()) {
            case AtomicShell.KFamily :
               res = KLMLine.KLMLineType.KEdge;
               break;
            case AtomicShell.LFamily :
               res = KLMLine.KLMLineType.LEdge;
               break;
            case AtomicShell.MFamily :
               res = KLMLine.KLMLineType.MEdge;
               break;
            case AtomicShell.NFamily :
               res = KLMLine.KLMLineType.NEdge;
               break;
            default :
               res = KLMLine.KLMLineType.InvalidType;
               break;
         }
         return res;
      }

      /**
       * compareTo - Allows to KLM lines to be ordered.
       *
       * @param kl
       *           KLMLine
       * @return int
       */
      @Override
      public int compareTo(KLMLine kl) {
         // Ordering: Transition, SumPeak, EscapePeak, Edge
         if (kl instanceof Edge edge) {
            return this.mShell.compareTo(edge.mShell);
         } else {
            // Edges are always last.
            return -1;
            // return (kl instanceof Transition) || (kl instanceof SumPeak) ||
            // (kl instanceof EscapePeak) || (kl instanceof Edge) ? -1 : 1;

         }
      }

      @Override
      public boolean equals(Object obj) {
         if (obj == null) {
            return false;
         }
         if (this == obj) {
            return true;
         }
         if (obj instanceof Edge edge) {
            return this.mShell == edge.mShell;
         }
         return false;
      }

   }

   public static class EscapePeak extends KLMLine implements Comparable<KLMLine> {
      private final XRayTransition mTransition;

      static private final XRayTransition SI_K = new XRayTransition(Element.Si, XRayTransition.KA1);

      public EscapePeak(XRayTransition xrt) throws EPQException {
         super(xrt.getEnergy() - SI_K.getEnergy(), 0.1);
         mTransition = xrt;
      }

      @Override
      public String toLabel(LabelType lt) {
         switch (lt) {
            case ELEMENT :
               return mTransition.getElement().toString() + " esc";
            case ELEMENT_ABBREV :
            case LARGE_ELEMENT :
               return mTransition.getElement().toAbbrev() + " esc";
            case SIEGBAHN :
               return mTransition.getSiegbahnName() + " esc";
            case IUPAC :
            case FAMILY :
               return mTransition.getElement().toAbbrev() + " esc";
            case NONE :
            default :
               return "";
         }
      }

      @Override
      public String toString() {
         return mTransition.toString() + " escape";
      }

      @Override
      public boolean contains(Element elm) {
         return mTransition.getElement().equals(elm);
      }

      @Override
      public AtomicShell getShell() {
         return mTransition.getDestination();
      }

      @Override
      public int hashCode() {
         return mTransition.hashCode() + 0xFEEE;
      }

      /**
       * Suggest a list of possible escape peaks for the specified element from
       * the specified list of possible lines with energy less then eMax
       *
       * @param elm
       * @return Set&lt;EscapePeak&gt;
       */
      static public Set<EscapePeak> suggestEscapePeak(Element elm) {
         int[] lines = {XRayTransition.KA1, XRayTransition.LA1, XRayTransition.MA1};
         final HashSet<EscapePeak> res = new HashSet<>();
         for (final int line : lines) {
            try {
               if (XRayTransition.exists(elm, line)) {
                  final XRayTransition xrt = new XRayTransition(elm, line);
                  final double e = xrt.getEnergy();
                  if (e > ToSI.keV(0.02) + SI_K.getEnergy()) {
                     res.add(new EscapePeak(xrt));
                  }
               }
            } catch (final Exception e) {
               // Just ignore it
            }
         }
         return res;
      }

      @Override
      public String toSiegbahn() {
         return mTransition.getSiegbahnName() + "-Si K";
      }

      @Override
      public KLMLineType getType() {
         return KLMLine.KLMLineType.EscapePeak;
      }

      /**
       * compareTo - Allows to KLM lines to be ordered.
       *
       * @param kl
       *           KLMLine
       * @return int
       */
      @Override
      public int compareTo(KLMLine kl) {
         // Ordering: Transition, SumPeak, EscapePeak, Edge
         if (kl instanceof EscapePeak esc) {
            return this.mTransition.compareTo(esc.mTransition);
         } else {
            return (kl instanceof Transition) || (kl instanceof SumPeak) ? -1 : 1;
         }
      }

      @Override
      public boolean equals(Object obj) {
         if (obj == null) {
            return false;
         }
         if (this == obj) {
            return true;
         }
         if (obj instanceof EscapePeak esc) {
            return this.mTransition == esc.mTransition;
         }
         return false;
      }

   }

   public static class SumPeak extends KLMLine implements Comparable<KLMLine> {
      private final ArrayList<XRayTransition> mTransitions;
      private final String mName;

      public SumPeak(XRayTransition xrt1, XRayTransition xrt2) throws EPQException {
         super(xrt1.getEnergy() + xrt2.getEnergy(), 0.1);
         mTransitions = new ArrayList<>();
         mTransitions.add(xrt1);
         mTransitions.add(xrt2);
         Collections.sort(mTransitions);
         if (xrt1.equals(xrt2)) {
            mName = "2\u00B7" + xrt1.toString();
         } else {
            mName = xrt1.toString() + "+" + xrt2.toString();
         }
      }

      public SumPeak(XRayTransition xrt1, XRayTransition xrt2, XRayTransition xrt3) throws EPQException {
         super(xrt1.getEnergy() + xrt2.getEnergy() + xrt3.getEnergy(), 0.1);
         mTransitions = new ArrayList<>();
         mTransitions.add(xrt1);
         mTransitions.add(xrt2);
         mTransitions.add(xrt3);
         Collections.sort(mTransitions);
         if (xrt1.equals(xrt2) && xrt1.equals(xrt3)) {
            mName = "3\u00B7" + xrt1.toString();
         } else {
            mName = xrt1.toString() + "+" + xrt2.toString() + "+" + xrt3.toString();
         }
      }

      private static double sumXRT(Collection<XRayTransition> xrts) {
         double res = 0.0;
         for (XRayTransition xrt : xrts) {
            try {
               res += xrt.getEnergy();
            } catch (EPQException e) {
               e.printStackTrace();
            }
         }
         return res;
      }

      public SumPeak(Collection<XRayTransition> xrts) throws EPQException {
         super(sumXRT(xrts), 0.1);
         TreeMap<XRayTransition, Integer> hm = new TreeMap<>();
         for (XRayTransition xrt : xrts) {
            hm.put(xrt, hm.getOrDefault(xrt, 0) + 1);
         }
         StringBuffer sb = new StringBuffer();
         for (Map.Entry<XRayTransition, Integer> me : hm.entrySet()) {
            if (!sb.isEmpty()) {
               sb.append("+");
            }
            if (me.getValue() > 1) {
               sb.append(me.getValue().toString());
               sb.append("\u00B7");
               sb.append(me.getKey().toString());
            } else {
               sb.append(me.getKey().toString());
            }
         }
         mTransitions = new ArrayList<>(xrts);
         Collections.sort(mTransitions);
         mName = sb.toString();
      }

      @Override
      public String toLabel(LabelType lt) {
         final StringBuffer sb = new StringBuffer();
         final TreeMap<XRayTransition, Integer> count = new TreeMap<>();
         for (final XRayTransition mTransition : mTransitions) {
            if (count.containsKey(mTransition)) {
               count.put(mTransition, Integer.valueOf(count.get(mTransition) + 1));
            } else {
               count.put(mTransition, Integer.valueOf(1));
            }
         }
         for (final Map.Entry<XRayTransition, Integer> me : count.entrySet()) {
            if (sb.length() > 0) {
               sb.append("+");
            }
            if (me.getValue().intValue() > 1) {
               sb.append(me.getValue().toString());
               sb.append("\u00B7");
            }
            final XRayTransition xrt = me.getKey();
            switch (lt) {
               case ELEMENT :
                  sb.append(xrt.getElement().toString());
                  break;
               case LARGE_ELEMENT :
               case ELEMENT_ABBREV :
                  sb.append(xrt.getElement().toAbbrev());
                  break;
               case SIEGBAHN :
                  sb.append(xrt.getSiegbahnName());
               case IUPAC :
                  sb.append(xrt.toString());
                  break;
               case FAMILY :
                  sb.append(xrt.getElement().toAbbrev() + " " + XRayTransitionSet.getClassNameForTransition(xrt.getTransitionIndex()));
                  break;
               case NONE :
               default :
                  continue;
            }
         }
         return sb.toString();
      }

      @Override
      public String toString() {
         return mName;
      }

      @Override
      public AtomicShell getShell() {
         return mTransitions.get(0).getDestination();
      }

      @Override
      public int hashCode() {
         int hash = 0x0;
         for (final XRayTransition xrt : mTransitions) {
            hash ^= xrt.hashCode();
         }
         return hash;
      }

      @Override
      public boolean contains(Element elm) {
         for (final XRayTransition xrt : mTransitions) {
            if (xrt.getElement().equals(elm)) {
               return true;
            }
         }
         return false;
      }

      static public Set<SumPeak> suggestSumPeaks(Set<Element> elms, double lowE, double highE, int order) {
         final HashSet<XRayTransition> xrts = new HashSet<>();
         for (final Element elm : elms) {
            for (XRayTransition xrt : new XRayTransitionSet(elm, 0.0, highE)) {
               if (xrt.getWeight(XRayTransition.NormalizeFamily) > 0.5) {
                  xrts.add(xrt);
               }
            }
         }
         final HashSet<SumPeak> res = new HashSet<>();
         for (XRayTransition xrt1 : xrts) {
            for (XRayTransition xrt2 : xrts) {
               try {
                  final double e1 = xrt1.getEnergy();
                  final double e2 = xrt2.getEnergy();
                  if (((e1 + e2) >= lowE) && ((e1 + e2) <= highE)) {
                     res.add(new SumPeak(xrt1, xrt2));
                  }
                  if (order > 2) {
                     for (XRayTransition xrt3 : xrts) {
                        final double e3 = xrt3.getEnergy();
                        if (((e1 + e2 + e3) >= lowE) && ((e1 + e2 + e3) <= highE)) {
                           res.add(new SumPeak(xrt1, xrt2, xrt3));
                        }
                     }
                  }
               } catch (final EPQException e) {
                  System.err.println("This should never happen!");
               }
            }
         }
         return res;
      }

      @Override
      public String toSiegbahn() {
         if (mTransitions.stream().allMatch(a -> a.equals(mTransitions.get(0)))) {
            return Integer.toString(mTransitions.size()) + "\u00B7" + mTransitions.get(0).getSiegbahnName();
         } else {
            return mTransitions.stream().map(xrt -> xrt.getSiegbahnName()).collect(Collectors.joining("+"));
         }
      }

      @Override
      public KLMLineType getType() {
         return KLMLine.KLMLineType.SumPeak;
      }

      /**
       * compareTo - Allows to KLM lines to be ordered.
       *
       * @param kl
       *           KLMLine
       * @return int
       */
      @Override
      public int compareTo(KLMLine kl) {
         // Ordering: Transition, SumPeak, EscapePeak, Edge
         if (kl instanceof SumPeak sum) {
            Iterator<XRayTransition> i1 = this.mTransitions.iterator();
            Iterator<XRayTransition> i2 = sum.mTransitions.iterator();
            while (i1.hasNext() && i2.hasNext()) {
               XRayTransition xrt1 = i1.next();
               XRayTransition xrt2 = i2.next();
               int cmp = xrt1.compareTo(xrt2);
               if (cmp != 0) {
                  return cmp;
               }
            }
            return i1.hasNext() ? -1 : (i2.hasNext() ? 1 : 0);
         } else {
            // return (kl instanceof Transition) || (kl instanceof SumPeak) ||
            // (kl instanceof EscapePeak) || (kl instanceof Edge) ? -1 : 1;
            return (kl instanceof Transition) ? -1 : 1;
         }
      }

      @Override
      public boolean equals(Object obj) {
         if (obj == null) {
            return false;
         }
         if (this == obj) {
            return true;
         }
         if (obj instanceof SumPeak sum) {
            Iterator<XRayTransition> i1 = mTransitions.iterator();
            Iterator<XRayTransition> i2 = sum.mTransitions.iterator();
            while (i1.hasNext() && i2.hasNext()) {
               if (i1.next() != i2.next()) {
                  return false;
               }
            }
            return !(i1.hasNext() || i2.hasNext());
         }
         return false;
      }
   }

   /**
    * Constructs a KLMLine representing a transition
    *
    * @param energy
    * @param amplitude
    */
   protected KLMLine(double energy, double amplitude) {
      mEnergy = energy;
      mAmplitude = amplitude;
   }

   /**
    * getEnergy - Returns the energy of a particular KLM line in Joules.
    *
    * @return double
    */
   public double getEnergy() {
      return mEnergy;
   }

   /**
    * getAmplitude - returns the amplitude of a particular KLM line.
    *
    * @return double
    */
   public double getAmplitude() {
      return mAmplitude;
   }

   abstract public boolean contains(Element elm);

   abstract public String toSiegbahn();

   public KLMLineType getType() {
      return KLMLine.KLMLineType.InvalidType;
   }

   static String exportKLMLines(Collection<KLMLine> lines) {
      StringBuffer sb = new StringBuffer(lines.size() * 100);
      for (KLMLine line : lines) {
         switch (line.getType()) {
            case InvalidType :
            case Satellite :
               // Skip
               continue;
            case SumPeak :
               sb.append("Sum:" + line.toString() + "\n");
               continue;
            case EscapePeak :
               sb.append("Escape:" + line.toString() + "\n");
               continue;
            case KEdge :
            case LEdge :
            case MEdge :
            case NEdge :
               sb.append("Edge:" + line.toString() + "\n");
               continue;
            case KTransition :
            case LTransition :
            case MTransition :
            case NTransition :
               sb.append("Transition:" + line.toString() + "\n");
               continue;
         }

      }
      return sb.toString();
   }

   static private KLMLine parseKLM(String str) throws EPQException {

      if (str.startsWith("Sum:")) {
         Collection<XRayTransition> lines = new ArrayList<>();
         String[] ss = str.substring(4).trim().split("\\+");
         for (String s : ss) {
            if (s.matches("[2-9]·.*")) {
               int n = Integer.parseInt(s.substring(0, 1));
               XRayTransition xrt = XRayTransition.parseString(s.substring(2).trim());
               for (int i = 0; i < n; ++i) {
                  lines.add(xrt);
               }
            } else {
               XRayTransition xrt = XRayTransition.parseString(s.trim());
               lines.add(xrt);
            }
         }
         return new KLMLine.SumPeak(lines);
      } else if (str.startsWith("Escape:")) {
         XRayTransition xrt = XRayTransition.parseString(str.substring(7, str.length() - 7).trim());
         return new KLMLine.EscapePeak(xrt);
      } else if (str.startsWith("Edge:")) {
         AtomicShell ass = AtomicShell.parseString(str.substring(5).trim());
         return new KLMLine.Edge(ass);
      } else if (str.startsWith("Transition:")) {
         XRayTransition xrt = XRayTransition.parseString(str.substring(11).trim());
         return new KLMLine.Transition(xrt);
      } else {
         throw new EPQException("Unrecognized type while parsing KLM lines.");
      }
   }

   static ArrayList<KLMLine> importKLMLines(Collection<String> lines) {
      ArrayList<KLMLine> res = new ArrayList<>();
      for (String line : lines) {
         try {
            if (!line.trim().isEmpty()) {
               res.add(KLMLine.parseKLM(line));
            }
         } catch (EPQException e) {
            e.printStackTrace();
         }
      }
      return res;
   }

}

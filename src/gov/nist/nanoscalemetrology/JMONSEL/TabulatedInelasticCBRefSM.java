/**
 *
 */
package gov.nist.nanoscalemetrology.JMONSEL;

import java.io.FileNotFoundException;
import java.util.Arrays;

import gov.nist.microanalysis.EPQLibrary.EPQFatalException;
import gov.nist.microanalysis.EPQLibrary.Material;
import gov.nist.microanalysis.NISTMonte.Electron;
import gov.nist.microanalysis.Utility.Math2;

/**
 * <p>
 * This is a temporary version of TabulatedInelasticSM. It contains some
 * differences that are under test. The reason for this second version is an
 * ambiguity in the meaning of PE kinetic energy in the DFT model. What is the
 * right kinetic energy to use? If E0 is the PE's kinetic energy when it is
 * outside the sample, then is the relevant kinetic energy for scattering = E0,
 * E0 - UCB, or E0 - USB, where UCB is the potential energy at the bottom of the
 * conduction band and UVB is the potential energy at the bottom of the
 * scattering band? The "scattering band" is the band that contains the
 * electrons from which the PE is scattering. I.e., in a conductor it is the
 * same as the conduction band. In an insulator it is the filled valence band,
 * the bottom of which is offset from the bottom of the conduction band by some
 * finite energy, generally something between 10 eV and 20 eV.
 * </p>
 * <p>
 * The ambiguity arises because the DFT model, which I presently use to compute
 * the scattering tables, is based ultimately on a dielectric function for a
 * free electron gas. Such a gas is a Fermi sea of electrons in which the least
 * energetic (bottom of the band) have 0 kinetic energy, i.e., the same energy
 * as a stationary primary electron. In reality, our electrons are always bound
 * within a crystal. I.e., they have a negative potential energy. The way this
 * is treated for conductors is simple. When the PE enters the crystal it gains
 * energy -UCB. It's new kinetic energy is referenced to the same 0 as the
 * electrons in the conduction band. That is, with this slightly higher kinetic
 * energy, we are back to the situation that was modeled by DFT, and it seems we
 * can appropriately apply it.
 * </p>
 * <p>
 * However, what if we are dealing with an insulator? In that case the target
 * electrons are in an even more tightly bound band, with band bottom a distance
 * Eoffset below the CB bottom. What should we do? Should we add this additional
 * amount to our estimate of PE kinetic energy (i.e., kinetic energy = E0 -
 * USB)? Does this restore the situation to that which was modeled by DFT? That
 * is the model implemented already in TabulatedInelasticSM.
 * </p>
 * <p>
 * Alternatively, should we continue to treat our PE as having kinetic energy E0
 * - UCB? That will be the model implemented here. A hand-waving argument for
 * this model is that the wave function of the PE can be written as a
 * superposition of eigenstates of the lattice (Bloch states) of the appropriate
 * energy. That energy is, however, high (higher even in most cases than
 * conduction band states). Low energy states have higher probability density
 * close to the atomic cores, where the potential energy is very low, but high
 * energy states, which must be orthogonal to those, are concentrated most
 * strongly in the regions between atomic cores, where the potential is flat and
 * close to UCB. Thus, the velocity characteristic of those states should be
 * approximately sqrt(2*(E0-UCB)/m). Even the bound valence electrons spend
 * considerable time in these regions. That's why the valence band is an
 * extended band rather than a localized atomic-core state.
 * </p>
 * <p>
 * On the other hand, here's a hand-waving argument for the former model: We are
 * interested in scattering events that transfer enough energy to overcome the
 * bandgap, and we are probably most interested in those that transfer enough to
 * free the SE entirely from the crystal. These are events that transfer
 * something on the order of 10 eV to 20 eV. Energy transfers of this size
 * suggest (in a classical picture) a small impact factor: PE and SE in close
 * proximity at the time of the energy transfer. Since the target electron has
 * characteristic potential energy USB, if our PE is very close to the target,
 * it should also have potential energy USB. Thus, even if on average the
 * electron is in parts of the sample characterized by potential energy UCB, at
 * the relevant time its potential energy is lower (and therefore its kinetic
 * energy higher).
 * </p>
 * <p>
 * I'm not sure which of these is the better model. This is the reason to
 * implement and try both.
 * </p>
 * <p>
 * The TabulatedInelasticSM differs from earlier implementations of
 * ScatterMechanism in that most of the properties of the scattering are not
 * hard coded into the object, but rather are dictated by input tables. So, for
 * example, the inverse mean free path for this mechanism is determined by
 * interpolating an IIMFP table associated with the material. The distribution
 * of energy losses and angle changes by the primary electron are similarly
 * determined by tables. In this way, almost all of the physics of the
 * scattering is contained in the tables, and the scatter mechanism object is
 * itself almost without physics content.
 * </p>
 * <p>
 * The tables are provided to the constructor via an array of strings, each of
 * which contains the full path to a file that can serve as input to construct a
 * NUTableInterpolation object. The required tables are documented in the
 * constructor notes. For 3 of the tables, the first input parameter is the
 * kinetic energy of the primary electron. By default, this kinetic energy is
 * assumed to be the PE's energy with respect to the bottom of the conduction
 * band. Sometimes, however, tables use a different convention. For example, in
 * a semiconductor or insulator tables may be computed for energies measured
 * with respect to the bottom of the valence band (the highest occupied band).
 * This can be accommodated by providing the constructor with the energy offset
 * between these two definitions of energy. That is, Eoffset = (energy of
 * conduction band bottom) - (the energy defined as the zero for purpose of the
 * tables).
 * </p>
 * <p>
 * There is, however, an exception to the above generalizations. There remains
 * some uncertainty in the literature over the connection between primary
 * electron energy loss/trajectory change on the one hand and secondary electron
 * (SE) final energy/final trajectory on the other. Part of the reason for this
 * difference lies in varying assumptions about the initial (pre-scattering)
 * energy and trajectory of the SE. Some of the different methods of treating
 * this connection have been implemented in this class. One of them must be
 * selected to instantiate an object of this class. The selection is made by
 * choosing a value for the methodSE parameter in the constructor. Following are
 * the methods that have been implemented:
 * </p>
 * <p>
 * methodSE = 1: This selection is an implementation of the method described by
 * Ding &amp; Shimizu in SCANNING 18 (1996) p. 92. If the PE energy loss, deltaE is
 * greater than a core level binding energy, the SE final energy is
 * deltaE-Ebinding. Otherwise, it is deltaE+EFermi, where EFermi is the Fermi
 * energy of the material. The final direction of the SE is determined from
 * conservation of momentum with the assumption that the SE initial momentum was
 * 0.
 * </p>
 * <p>
 * methodSE = 2: This selection is an implementation of the method described by
 * Ding, Tang, &amp; Shimizu in J.Appl.Phys. 89 (2001) p. 718. If deltaE is greater
 * than a core level binding energy the treatment is the same as methodSE = 1.
 * If not, the SE final energy is deltaE + E'. If E' were the Fermi energy this
 * would be the same as methodSE = 1. However, E' lies in the range max(0,EFermi
 * - deltaE) &le; E' &le; EFermi. The value of E' is determined probabilistically
 * based upon the free electron densities of occupied and unoccupied states.
 * </p>
 * <p>
 * methodSE = 3: This selection is my modified version of the method described
 * by Mao et al. in J.Appl.Phys. 104 (2008) article #114907. The scattering
 * event is assigned as either an electron-electron event or an electron-plasmon
 * event based upon the momentum transfer. Plasmon events are treated as in
 * methodSE = 2, except that the SE direction is isotropic. Electron-electron
 * events have energy and direction determined as described by Mao et al.
 * </p>
 * <p>
 * Interpolation is implemented using the NUTableInterpolation class. That class
 * natively allows extrapolation, but TabulatedInelasticSM checks input
 * parameters and forbids it. This means the user-provided tables must cover the
 * full range of energies that will be encountered in the simulation. This
 * generally means the tables that take PE energy as one of the inputs should
 * cover energies from the material's Fermi energy up to at least the energy of
 * electrons generated by the electron gun.
 * </p>
 * <p>
 * This scattering mechanism uses the following material properties, which
 * therefore need to be properly defined: coreEnergy array, energyCBbottom,
 * bandgap, and workfunction.
 * </p>
 * <p>
 * Copyright: Pursuant to title 17 Section 105 of the United States Code this
 * software is not subject to copyright protection and is in the public domain.
 * </p>
 * <p>
 * Company: National Institute of Standards and Technology
 * </p>
 *
 * @author John Villarrubia
 * @version 1.0
 */

public class TabulatedInelasticCBRefSM
   extends
   ScatterMechanism {

   private final int methodSE;
   private double energyOffset = 0.;

   private NUTableInterpolation tableIIMFP;
   private NUTableInterpolation tableReducedDeltaE;
   private NUTableInterpolation tableTheta;
   private NUTableInterpolation tableSEE0;
   /*
    * The Fermi energy is defined for SEmaterials as the position of the highest
    * occupied state relative to the bottom of the conduction band. This makes
    * it equal to what we usually think of as the Fermi energy for a metal and
    * equal to -bandgap for an insulator or semiconductor. In my scattering
    * table notes, on the other hand, Fermi energy was the distance from bottom
    * to highest occupied state in whatever band was responsible for scattering.
    * This makes no change for metals, but it changes the value when bandgap !=
    * 0. offsetFermiEnergy (next line) corresponds to my Mathematica notebook
    * definition.
    */
   private double offsetFermiEnergy;
   private double energyCBbottom;
   private double minEgenSE = 0.;
   private double workfunction;
   private double bandgap;
   private double energyGap;
   private boolean defaultRatios = true;
   private double[][] cumulativeBranchingProbabilities = null;

   /*
    * bEref is the energy (relative to conduction band bottom) to which core
    * level binding energies are referenced. This is generally the Fermi energy
    * for metals and 0 for insulators or semiconductors.
    */
   private double bEref;
   private final double[] kEa = new double[1]; // For convenience, because 1-d
   // tables
   // still require an array for input
   private Double[] coreEnergies;
   private final double[] interpInput = new double[3];

   // Allowed energy ranges for interpolation table inputs
   private double[] tableEiDomain;
   private double[] tableIIMFPEiDomain;
   // Range of allowed SE initial energies on output
   private double[] energyRangeSE0;

   // temp? IIMFP multiplier
   private double rateMult = 1.;

   /*
    * E0 method. E0 is the energy available to ionize an inner shell. My
    * original method based on the description in Ding &amp; Shimizu's SCANNING
    * article was to assume that deltaE, i.e., all energy transferred in an
    * inelastic event, is available for ionization. Ding et al later modified
    * this procedure. By the plasmon dispersion (I use Penn's form) deltaE can
    * be broken down into a q=0 part and a q != 0 part. The q = 0 part is E0.
    * When E0fromDispersion = false I used deltaE. When it is true I use Ding's
    * later method.
    */
   private boolean E0fromDispersion = false;

   /**
    * Constructs a TabulatedInelasticCBRefSM for the specified material.
    *
    * @param mat - a SEmaterial that is the material within which scattering
    *           occurs.
    * @param methodSE - an int that determines the method and assumptions by
    *           which SE energies and angles are determined in a scattering
    *           event. See the description in the class documentation.
    * @param tables - an array of strings. The strings contain the full paths
    *           and file names of the interpolation tables in this order:
    *           String[0] = the IIMFP table (inverse inelastic mean free path
    *           vs. primary electron energy (EO), String[1] = the reduced deltaE
    *           table (deltaE/E0 vs E0 and r, with r a random number), table[2]
    *           = the theta table (scattering angle of PE vs. E0,
    *           deltaE/(E0-EFermi), r) and table[3] = the table of SE initial
    *           energy vs. deltaE and r.
    */
   public TabulatedInelasticCBRefSM(SEmaterial mat, int methodSE, String[] tables) {
      this(mat, methodSE, tables, 0.);
   }

   /**
    * <p>
    * Constructs a TabulatedInelasticCBRefSM for the specified material. This
    * form of the constructor has an additional argument, energyOffset, allowing
    * this parameter to be set to a value other than its default value of 0.
    * This parameter is still needed in the present model, despite that this
    * model does not use it to estimate the PE's kinetic energy, because Eoffset
    * affects the SE kinetic energy after energy transfer.
    * </p>
    * <p>
    * energyOffset = (energy of conduction band bottom) - (the energy defined as
    * the zero for purpose of the tables, generally the scattering band bottom)
    */
   public TabulatedInelasticCBRefSM(SEmaterial mat, int methodSE, String[] tables, double energyOffset) {
      super();
      if((methodSE != 2) && (methodSE != 3))
         methodSE = 1; // Make sure methodSE is valid
      this.methodSE = methodSE;

      /* Read interpolation tables into memory */
      try {
         tableIIMFP = NUTableInterpolation.getInstance(tables[0]);
      }
      catch(final FileNotFoundException e1) {
         throw new EPQFatalException("File " + tables[0] + " not found.");
      }
      try {
         tableReducedDeltaE = NUTableInterpolation.getInstance(tables[1]);
      }
      catch(final FileNotFoundException e1) {
         throw new EPQFatalException("File " + tables[1] + " not found.");
      }
      try {
         tableTheta = NUTableInterpolation.getInstance(tables[2]);
      }
      catch(final FileNotFoundException e1) {
         throw new EPQFatalException("File " + tables[2] + " not found.");
      }
      if((methodSE == 2) || (methodSE == 3))
         try {
            tableSEE0 = NUTableInterpolation.getInstance(tables[3]);
            energyRangeSE0 = tableSEE0.getRange();
         }
         catch(final FileNotFoundException e1) {
            throw new EPQFatalException("File " + tables[3] + " not found.");
         }

      this.energyOffset = energyOffset;
      setMaterial(mat);
   }

   /*
    * (non-Javadoc)
    * @see
    * gov.nist.nanoscalemetrology.JMONSEL.ScatterMechanism#scatter(gov.nist.
    * microanalysis.NISTMonte.Electron)
    */
   @Override
   public Electron scatter(Electron pe) {

      final double kE = pe.getEnergy(); // PE initial energy rel to CB bottom

      if(kE < tableEiDomain[0])
         /*
          * This might happen if something, e.g., electrostatic potential
          * difference, reduces electron energy between the time we determine
          * that scattering happens (at the beginning of a step) and the time it
          * actually occurs (at the end of the step). In this case, we simply
          * don't scatter.
          */
         return null;
      if(kE > tableEiDomain[1])
         throw new EPQFatalException("PE energy " + Double.toString(kE) + " is outside the interpolation table interval of ["
               + Double.toString(tableEiDomain[0]) + "," + Double.toString(tableEiDomain[1]) + "]");
      double theta = 0.;
      double phi = 0.; // PE trajectory parameters
      double energySE, thetaSE, phiSE; // SE trajectory parameters

      final double[] randoms = new double[] {
         Math2.rgen.nextDouble(),
         Math2.rgen.nextDouble(),
         Math2.rgen.nextDouble(),
         Math2.rgen.nextDouble()
      };
      interpInput[0] = kE;
      interpInput[1] = randoms[0];
      // Energy loss by PE
      double deltaE = kE * tableReducedDeltaE.interpolate(interpInput, 3);
      /*
       * Cubic interpolation of the table can undershoot. Treat deltaE close to
       * but below the energyGap as such undershoot and correct it.
       */
      if((deltaE < energyGap) && (deltaE > (0.95 * energyGap)))
         deltaE = energyGap;
      /*
       * Larger discrepancies are most likely because we've been supplied an
       * empirical table that includes non-electronic energy losses (e.g.,
       * scattering from phonons). These should really be handled separately
       * because our model is only valid for electrons & plasmons. (E.g., phonon
       * of energy deltaE carries very different momentum from electron of
       * energy deltaE, so scattering angles can't be determined in the present
       * model.) We skip the angular scattering part for such events. Any
       * generated SE will be in the bandgap, so most likely dropped anyway. We
       * return after we deal with the PE energy loss.
       */

      final double theta0PE = pe.getTheta(); // Remember original direction;
      final double phi0PE = pe.getPhi(); // to use for SE
      if(deltaE >= bandgap) {
         // Determine theta and phi here
         /*
          * First, the reduced energy. This parameter ranges from 0 to 1 as
          * deltaE ranges from its minimium to maximum value.
          */
         interpInput[1] = (deltaE - energyGap) / (kE - energyGap);
         /*
          * The reduced energy can on rare occasions, as a result of
          * interpolation error, lie slightly outside its physically determined
          * interval of [0,1]. If it does, clip it to the boundary.
          */
         if(interpInput[1] > 1.)
            interpInput[1] = 1.;
         else if(interpInput[1] < 0.)
            interpInput[1] = 0.;
         interpInput[2] = randoms[1];
         theta = tableTheta.interpolate(interpInput, 3);
         phi = 2. * Math.PI * randoms[2];
         /*
          * Update PE trajectory. Note that the energy of the PE is decremented
          * by deltaE. Any continuous energy loss formula should only account
          * for losses not including this one.
          */
         pe.updateDirection(theta, phi);
      }

      pe.setEnergy(kE - deltaE);

      // Determine SE final energy and trajectory
      Electron se = null;
      double be = 0.;

      /*
       * Some measured ELF data have nonzero values for deltaE less than the
       * bandgap, and it is consequently possible that some scattering tables
       * that retain this part of the data will include such loss events. These
       * may, for example, correspond to phonon losses. They presumably do not
       * correspond to generation of mobile SE, since there are no empty mobile
       * states in the gap. We therefore return no SE for such events.
       */
      if(deltaE < bandgap)
         return null;

      final double Eq = (2. * kE) - deltaE - (2. * Math.sqrt(kE * (kE - deltaE)) * Math.cos(theta));

      switch(methodSE) {
         case 1:
            /*
             * In the following formula, offsetFermiEnergy - energyOffset is the
             * Fermi energy re-referenced to the bottom of the conduction band.
             * If b (the nearest lower core level binding energy) is zero, then
             * this is the electron's initial energy. (mode 1 assumes target
             * electrons come from the Fermi level.) If b>0 then the Fermi
             * energy - b is still the electron's initial energy, since the core
             * level energies are referenced to the Fermi level. Either way,
             * adding deltaE gives the SE's final energy.
             */
            energySE = (deltaE + bEref) - pickBE(Eq, deltaE);
            if((energySE + energyCBbottom) < minEgenSE)
               return null;
            thetaSE = (Math.PI / 2.) - theta;
            phiSE = phi + Math.PI;
            // Generate SE, apply energy loss and trajectory change to SE here
            se = new Electron(pe, theta0PE, phi0PE, energySE);
            se.updateDirection(thetaSE, phiSE);
            break;
         case 2:
            be = pickBE(Eq, deltaE);
            if(be > 0.)
               energySE = (deltaE + bEref) - be;
            else {
               interpInput[0] = deltaE;
               interpInput[1] = randoms[3];
               double energy0SE = tableSEE0.interpolate(interpInput, 3);
               /*
                * The values in the SEE0 table should range from 0 to EFermi,
                * which represents the range of allowed values. If the
                * interpolated value overshoots, clip it.
                */
               if(energy0SE < energyRangeSE0[0])
                  energy0SE = energyRangeSE0[0];
               else if(energy0SE > energyRangeSE0[1])
                  energy0SE = energyRangeSE0[1];
               energySE = (deltaE + energy0SE) - energyOffset;
            }
            if((energySE + energyCBbottom) < minEgenSE)
               return null;
            thetaSE = (Math.PI / 2.) - theta;
            phiSE = phi + Math.PI;
            // Generate SE, apply energy loss and trajectory change to SE here
            se = new Electron(pe, theta0PE, phi0PE, energySE);
            se.updateDirection(thetaSE, phiSE);
            break;
         case 3:
            be = pickBE(Eq, deltaE);
            if(be > 0.) { // core level excitation
               energySE = (deltaE + bEref) - be;
               if((energySE + energyCBbottom) < minEgenSE)
                  return null;
               /*
                * I'm going to approximate the angle distribution as isotropic
                * for now.
                */
               thetaSE = Math.acos(1. - (2. * Math2.rgen.nextDouble()));
               phiSE = 2. * Math.PI * Math2.rgen.nextDouble();
               // Generate SE, apply energy loss and trajectory change to SE
               // here
               se = new Electron(pe, thetaSE, phiSE, energySE);
            } else { // SE generation from extended band
               final double root = 2. * Math.sqrt(offsetFermiEnergy * (offsetFermiEnergy + deltaE));
               final double sum = (2. * offsetFermiEnergy) + deltaE;
               final double Eqmin = sum - root;
               final double Eqmax = sum + root;
               if((Eqmin <= Eq) && (Eq <= Eqmax)) { // single-electron
                  // scattering
                  final double[] energytheta = simESEf(Eq, deltaE, randoms[3]);
                  energySE = energytheta[0] - energyOffset;
                  if((energySE + energyCBbottom) < minEgenSE)
                     return null;
                  // Generate SE in PE direction with correct energy
                  se = new Electron(pe, theta0PE, phi0PE, energySE);

                  // Determine angles of SE q vector relative to PE original
                  // direction
                  thetaSE = (Math.PI / 2.) - theta;
                  phiSE = phi + Math.PI;

                  // Combine with adjustment for additional simESEf deflection

                  final double[] newdir = updateDirection(thetaSE, phiSE, energytheta[1], 2. * Math.PI
                        * Math2.rgen.nextDouble());
                  // Update SE direction by this combined amount
                  se.updateDirection(newdir[0], newdir[1]);

               } else { // plasmon scattering
                  interpInput[0] = deltaE;
                  interpInput[1] = randoms[3];
                  double energy0SE = tableSEE0.interpolate(interpInput, 3);
                  /*
                   * The values in the SEE0 table should range from 0 to EFermi,
                   * which represents the range of allowed values. If the
                   * interpolated value overshoots, clip it.
                   */
                  if(energy0SE < energyRangeSE0[0])
                     energy0SE = energyRangeSE0[0];
                  else if(energy0SE > energyRangeSE0[1])
                     energy0SE = energyRangeSE0[1];
                  energySE = (deltaE + energy0SE) - energyOffset;
                  if((energySE + energyCBbottom) < minEgenSE)
                     return null;
                  /*
                   * For plasmon scattering, mode 3 assumes the plasmon
                   * "forgets" the momentum of the event that created it before
                   * it decays into an electron-hole pair. The angular
                   * distribution is therefore isotropic.
                   */
                  thetaSE = Math.acos(1. - (2. * Math2.rgen.nextDouble()));
                  phiSE = 2 * Math.PI * Math2.rgen.nextDouble();
                  // Generate SE, apply energy loss and trajectory change to SE
                  // here
                  se = new Electron(pe, thetaSE, phiSE, energySE);
               }
            }
            break;
         default:
            se = null;
            break;
      }

      return se;
   }

   /**
    * Updates a direction theta, phi by dtheta, dphi. This is the same algorithm
    * used by the Electron class to deflect an electron, except that it accepts
    * the initial angles in addition to the deflections and it returns the final
    * angles.
    *
    * @param theta double - The original polar angle
    * @param phi double - The original azimuthal angle
    * @param dTheta double - The deflection polar angle (0 = no deflection)
    * @param dPhi double - The deflection azimuthal angle
    */
   private double[] updateDirection(double theta, double phi, double dTheta, double dPhi) {

      final double ct = Math.cos(theta), st = Math.sin(theta);
      final double cp = Math.cos(phi), sp = Math.sin(phi);
      final double ca = Math.cos(dTheta), sa = Math.sin(dTheta);
      final double cb = Math.cos(dPhi);

      final double xx = (cb * ct * sa) + (ca * st);
      final double yy = sa * Math.sin(dPhi);
      final double dx = (cp * xx) - (sp * yy);
      final double dy = (cp * yy) + (sp * xx);
      final double dz = (ca * ct) - (cb * sa * st);

      theta = Math.atan2(Math.sqrt((dx * dx) + (dy * dy)), dz);
      phi = Math.atan2(dy, dx);
      return new double[] {
         theta,
         phi
      };
   }

   /*
    * simESEf is a private utility that computes the final SE energy for single
    * electron collisions. It also returns the polar angle of the final SE
    * trajectory in a reference frame with z axis in the direction of momentum
    * transfer, q. The inputs are Eq (the energy obtained from (hbar q)^2/(2m)
    * where q is the momentum transfer in the single electron event, deltaE is
    * the energy transferred in the event, and r is a random number from 0 to 1.
    * My derivation (in SimulatingALADingShimizu.nb) follows the one in Mao et
    * al., J. Appl. Phys. 2009.
    */
   private double[] simESEf(double Eq, double deltaE, double r) {
      final double q = Math.sqrt(Eq);
      final double kz = (deltaE - Eq) / 2. / q;
      final double kzf = kz + q;
      final double Ezq = kzf * kzf;
      double minE = (offsetFermiEnergy + bandgap) - Ezq;
      if(minE < 0.)
         minE = 0.;
      final double maxE = offsetFermiEnergy - (kz * kz);
      assert minE <= maxE;
      final double Exy = (minE * (1. - r)) + (maxE * r);
      final double ESEf = Exy + Ezq;
      final double theta = Math.acos(kzf / Math.sqrt(ESEf));
      return new double[] {
         ESEf,
         theta
      };
   }

   /*
    * This is a private utility used to determine the binding energy associated
    * with the secondary electron excitation channel. The default method
    * inherited from Ding &amp; Shimizu's SCANNING paper is to simply chose the
    * highest binding energy that is lower than deltaE. If the user has
    * specified branching ratios to associate with the various binding energies,
    * the alternative method chooses a binding energy at random with probability
    * consistent with the supplied ratios.
    */

   private double pickBE(double Eq, double deltaE) {
      int i;
      /*
       * Detect and return immediately in the most common case (deltaE too small
       * for inner shell excitation)
       */
      if((coreEnergies.length > 0) && (deltaE <= coreEnergies[0]))
         return 0.;
      /*
       * Arrive here if there is enough energy in principle to free an inner
       * shell electron. In this case we must compute E0 from the dispersion
       * equation, given Eq and deltaE.
       */
      double energyForIonization;
      if(E0fromDispersion)
         energyForIonization = computeE0fromDispersion(Eq, deltaE);
      else
         energyForIonization = deltaE;

      for(i = 0; (i < coreEnergies.length) && (coreEnergies[i] <= energyForIonization); i++)
         ;
      if(i == 0)
         return 0.;
      if(defaultRatios)
         /*
          * The advertised default behavior, as described by Ding & Shimizu
          * (Scanning).
          */
         return coreEnergies[i - 1];
      else {
         final double[] cprob = cumulativeBranchingProbabilities[i - 1];
         final double r = Math2.rgen.nextDouble();
         int index = Arrays.binarySearch(cprob, r);
         /*
          * Above binary search returns a positive index in the rare case when r
          * is exactly in the table. In such a case the energy we want is
          * coreEnergies[index]. Usually, though, there is no exact match for r.
          * In such case binarySearch returns index=-insertionPoint-1, so
          * insertionPoint (the index of the first table value that is larger
          * than r) is -(index+1). The energy we want though, is the next
          * smaller one than this, i.e., the one at -(index+1)-1 = -index-2. If
          * this remains less than 0, it means ALL table entries were greater
          * than r, in which case the energy we want is 0.
          */

         if(index < 0) {
            index = -2 - index;
            if(index < 0)
               return 0.;
         }

         return coreEnergies[index];
      }
   }

   /**
    * @param Eq
    * @param deltaE
    * @return
    */
   private double computeE0fromDispersion(double Eq, double deltaE) {
      /*
       * Note to self: The derivation of this algorithm is in
       * Y:\proj\linewidth\jvillar
       * \develop\NewMONSEL\Physics\DielectricDevelopment\
       * SimulatingALADingShimizu.nb
       */
      if(Eq == 0.)
         return deltaE;
      /* Precompute quantities we're going to need more than once. */
      final double x = Eq / deltaE;
      /*
       * The numeric constant on the next line is a constant that appears in the
       * solution of Penn's dispersion equation. It has units of 1/Joule.
       */
      final double y = 3.77300614251479e17 * Eq;
      final double x13 = Math.pow(x, 1 / 3.);
      final double x2 = x * x;
      final double c1 = x2 * x2 * (27. + (18. * y) + (2 * y * y));
      final double c2 = x2 * (27. + (4. * y));
      final double c3 = c2 - 27.;
      final double c4 = x2 * (3. + y);
      final double c6 = 1 - x2;
      if(c3 > 0.) {
         final double tanPart = Math.atan2(3. * c6 * Math.sqrt(3. * c3 * c6), (18. * c4) - c1 - 27.);
         double trigPart = Math.cos(tanPart / 3.);
         final double prefactor = 2. * x * Math.sqrt(y * (y - (c6 * (y + 6.))));
         trigPart *= prefactor;
         final double result = deltaE * Math.sqrt(((3. - c4) + trigPart) / 3.);
         return result;
      } else {
         final double c5 = Math.pow(-c1 + (18. * c4) + (3. * (-9. + (c6 * Math.sqrt(-(3. * c3 * c6))))), 1. / 3.);
         final double term1 = -2. * c4;
         final double x23 = x13 * x13;
         final double x43 = x * x13;
         final double y13 = Math.pow(y, 1. / 3.);
         final double y23 = y13 * y13;
         final double x43y23 = x43 * y23;
         final double two13 = Math.pow(2., 1. / 3.);
         final double term2 = -12. * two13 * x43y23 * c6;
         final double term3 = 2. * two13 * x43y23 * x2 * y;
         final double term23 = (term2 + term3) / c5;
         final double term4 = two13 * two13 * x23 * y13 * c5;
         final double result = deltaE * Math.sqrt((6 + term1 + term23 + term4) / 6.);
         return result;
      }
   }

   /*
    * (non-Javadoc)
    * @see
    * gov.nist.nanoscalemetrology.JMONSEL.ScatterMechanism#scatterRate(gov.nist
    * .microanalysis.NISTMonte.Electron)
    */
   @Override
   public double scatterRate(Electron pe) {
      kEa[0] = pe.getEnergy(); // The PE kinetic energy
      /*
       * The PE kinetic energy can fall below the minimum in the table for
       * materials with a energyGap. In this case the actual scatter rate is 0.
       */
      if(kEa[0] < tableIIMFPEiDomain[0])
         return 0.;
      if(kEa[0] > tableIIMFPEiDomain[1])
         throw new EPQFatalException("PE energy " + Double.toString(kEa[0]) + " exceeds interpolation table maximum energy of "
               + Double.toString(tableIIMFPEiDomain[1]));

      /*
       * I do only first order interpolation below because I noticed for some
       * tables that I get negative interpolated values. This happens despite
       * having all positive values in the table, because the scatter rate
       * approaches 0 at the Fermi level, leaving open the possibility of
       * overshoot. Possible approaches to avoid overshoot are to use linear
       * interpolation or to clip the result (as I do below for other tables).
       * Clipping to 0 seems a bad choice here, because it results in an
       * infinite inelastic free path.
       */
      final double result = rateMult * tableIIMFP.interpolate(kEa, 1);
      return result;
   }

   /*
    * (non-Javadoc)
    * @see
    * gov.nist.nanoscalemetrology.JMONSEL.ScatterMechanism#setMaterial(gov.nist
    * .microanalysis.EPQLibrary.Material)
    */
   @Override
   public void setMaterial(Material mat) {
      if(!(mat instanceof SEmaterial))
         throw new EPQFatalException("Material " + mat.toString()
               + " is not an SEmaterial as required for TabulatedInelasticSM.");

      final SEmaterial semat = (SEmaterial) mat;

      energyCBbottom = semat.getEnergyCBbottom();
      workfunction = semat.getWorkfunction();
      bandgap = semat.getBandgap();
      energyGap = bandgap;
      offsetFermiEnergy = semat.getEFermi() + energyOffset;
      coreEnergies = semat.getCoreEnergyArray();
      /*
       * My source for binding energies provides them relative to vacuum for a
       * few nobel and binary gases that we're unlikely to use, relative to the
       * Fermi level for metals, and relative to the top of the valence band for
       * insulators and semiconductors. bEref gives us our offsets to bottom of
       * conduction band.
       */
      if(bandgap > 0.)
         bEref = -bandgap;
      else
         bEref = semat.getEFermi();

      /*
       * tableEiDomain must be the energy range that is valid for *all* required
       * tables that take the PE initial energy as an input parameter.
       */
      tableEiDomain = tableReducedDeltaE.getDomain()[0];
      final double[] thetaTableEiDomain = tableTheta.getDomain()[0];
      if(thetaTableEiDomain[0] > tableEiDomain[0])
         tableEiDomain[0] = thetaTableEiDomain[0];
      if(thetaTableEiDomain[1] < tableEiDomain[1])
         tableEiDomain[1] = thetaTableEiDomain[1];
      tableIIMFPEiDomain = tableIIMFP.getDomain()[0];
   }

   /**
    * setMinEgenSE -- Sets the minimum energy for generated SE. Default is 0.
    * This model measures this energy relative to vacuum = 0. That is, minEgenSE
    * represents the energy that the SE would have if it immediately escaped the
    * sample (overcoming the work function barrier) without further loss of
    * energy. The default value of 0 means the lowest energy SE that are
    * generated are those that have barely enough energy to escape. It can be
    * set higher than this to turn off SE with higher energies, if they are not
    * relevant for a particular simulation. Note that the minEgenSE definition
    * here differs from the MollerInelasticSM model, which follows the
    * definition in the original MONSEL series. According to that definition,
    * the minEgenSE referred to the SE energy inside the sample. Accordingly, it
    * differed from this definition by an amount equal to the work function.
    *
    * @param minEgenSE The minEgenSE to set.
    */
   public void setMinEgenSE(double minEgenSE) {
      if(minEgenSE > -workfunction)
         this.minEgenSE = minEgenSE;
      else
         throw new EPQFatalException("Illegal minEgenSE.");
   }

   /**
    * @return Returns the minEgenSE.
    */
   public double getMinEgenSE() {
      return minEgenSE;
   }

   /**
    * Returns the value of the SE method parameter that is being used. This is
    * normally the same as the value supplied to the constructor.
    *
    * @return the methodSE
    */
   public int getMethodSE() {
      return methodSE;
   }

   /**
    * <p>
    * Branching ratios control how this class associates a core (binding) energy
    * with an excitation. If deltaE is the energy lost by the primary electron
    * in a scattering event, the secondary electron's final energy is equal to
    * its initial energy plus deltaE, but what is its initial energy? Was it
    * originally in the valence band, or in a deeper bound state? Only initial
    * states with binding energies less than deltaE contribute a channel for
    * energy loss. Typically the energy loss function (ELF) has a discontinuous
    * increase when deltE is equal to its energy. TabulatedInelasticSM assumes
    * the excitation channel that does not involve one of these bound electrons
    * is unaffected by the presence of the new excitation channel. Thus, the
    * ratio, r, of the ELF value from the left to the ELF value from the right
    * represents the probability that the new channel is not involved, and 1-r
    * is the probability that it is. These r ratios can easily be determined
    * from the ELF vs. energy curve. ratios (the argument of this method) is an
    * array of these ratios. The first ratio in the array is associated with the
    * lowest nonzero core energy (i.e., the first non-valence band bound state).
    * The length of the array must be equal to the length of the material's
    * coreEnergies array.
    * </p>
    * <p>
    * The default behavior, if this method is not called or if it is called with
    * no argument, is to assume all entries are 0. That is, the largest eligible
    * binding energy state is assumed to be the one associated with the
    * excitation channel. This is the method described by Ding &amp;Shimizu in
    * SCANNING.
    * </p>
    */
   public void setBranchingRatios() {
      defaultRatios = true;
   }

   public void setBranchingRatios(double[] ratios) {
      defaultRatios = false;
      final int cElen = coreEnergies.length;
      if(ratios.length != cElen)
         throw new EPQFatalException("The number of branching ratios must be equal to the number of core energies, "
               + Integer.toString(cElen) + " in this case.");
      final double[][] probabilities = new double[cElen][];
      cumulativeBranchingProbabilities = new double[cElen][];
      probabilities[0] = new double[] {
         ratios[0]
      };
      cumulativeBranchingProbabilities[0] = new double[] {
         ratios[0]
      };
      for(int i = 1; i < cElen; i++) {
         probabilities[i] = new double[i + 1];
         cumulativeBranchingProbabilities[i] = new double[i + 1];
         for(int j = 0; j < i; j++)
            probabilities[i][j] = probabilities[i - 1][j] * ratios[i];
         probabilities[i][i] = (1. - ratios[i - 1]) * ratios[i];
         cumulativeBranchingProbabilities[i][0] = probabilities[i][0];
         for(int j = 1; j <= i; j++)
            cumulativeBranchingProbabilities[i][j] = cumulativeBranchingProbabilities[i][j - 1] + probabilities[i][j];
      }
   }

   /**
    * <p>
    * This method was added to deal with LiF and similar materials. The
    * distinction between energyGap and bandgap is this: JMONSEL understands the
    * bandgap to be the distance between the top of the valence band and the
    * bottom of the conduction band. The energyGap is the value of the smallest
    * allowed deltaE in the scattering tables.
    * </p>
    * <p>
    * These two ordinarily are the same, and they are set equal by default.
    * However, they can differ if there are significant bound states within the
    * gap, scattering into which is included in the scattering tables. In LiF,
    * for instance, there are excitons just below the bottom of the conduction
    * band. These were included in the energy loss function from which the
    * scattering tables were derived in order to include their effect on the
    * stopping power.
    * </p>
    * <p>
    * The rule is, energyGap should be equal to the value used to compute the
    * scattering tables, while bandgap, which is used to determine the position
    * of the 0 of kinetic energy for free electrons in the crystal, should be
    * the distance between the bands. If these are the same (usual case) this
    * method need not be called. If they differ, use this method to distinguish
    * them.
    * </p>
    *
    * @param energyGap
    * @return
    */
   public void setEnergyGap(double energyGap) {
      this.energyGap = energyGap;
   }

   /**
    * The scatterRate or IIMFP (inverse inelastic mean free path) derived from
    * the associated table is multiplied by the factor provided via this routine
    * (default = 1). I'm not sure if I'm going to keep this. I added it as an ad
    * hoc way to simulate the effect of porosity in a material. The multiplier
    * can be set equal to the fill fraction, i.e., the fraction of the space
    * actually occupied by the material, in order to simulate the effect of many
    * tiny pores in the material.
    *
    * @param rateMult
    */
   public void setRateMult(double rateMult) {
      this.rateMult = rateMult;
   }

   /**
    * Gets the current value assigned to E0fromDispersion
    *
    * @return Returns the E0fromDispersion.
    */
   public boolean isE0fromDispersion() {
      return E0fromDispersion;
   }

   /**
    * Sets the value assigned to E0fromDispersion. If E0fromDispersion = false
    * (the default) JMONSEL continutes to use its original method (Ding &
    * Shimizu's SCANNING method) for determining the energy available to ionize
    * an inner shell. This method assumes the shell may be ionized whenever
    * deltaE (the energy transferred to the SE in an inelastic event) is greater
    * than the ionization energy. If E0fromDispersion = true, it uses Ding et
    * al's later method, in which 0-momentum part (E0) of the energy is computed
    * from the plasmon dispersion for an event which transfers deltaE. Inner
    * shell ionization can only happen if E0 &gt; ionization energy. This is more
    * restrictive than deltaE &gt; ionization energy.
    *
    * @param e0fromDispersion The value to which to set E0fromDispersion.
    */
   public void setE0fromDispersion(boolean e0fromDispersion) {
      E0fromDispersion = e0fromDispersion;
   }

}

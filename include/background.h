/** @file background.h Documented includes for background module */

#ifndef __BACKGROUND__
#define __BACKGROUND__

#include "common.h"
#include "quadrature.h"
#include "growTable.h"
#include "arrays.h"
#include "dei_rkck.h"
#include "parser.h"


// #include "complex.h" //kept here in case we need complex numbers one day
//we take this gamma function from gsl.
// static double complex gsl_sf_gamma(double complex z) {
static double gsl_sf_gamma(double z) {
    /* Lanczos coefficients for g = 7 */
    static double p[] = {
        0.99999999999980993227684700473478,
        676.520368121885098567009190444019,
       -1259.13921672240287047156078755283,
        771.3234287776530788486528258894,
       -176.61502916214059906584551354,
        12.507343278686904814458936853,
       -0.13857109526572011689554707,
        9.984369578019570859563e-6,
        1.50563273514931155834e-7
    };

    // if(creal(z) < 0.5)
    //     return M_PI / (sin(M_PI*z)*gsl_sf_gamma(1. - z));

    z -= 1;
    // double complex x = p[0];
    double x = p[0];
    for(int n = 1; n < 9; n++)
        x += p[n] / (z + 1.*n);
    // double complex t = z + 7.5;
    double t = z + 7.5;
    return sqrt(2*M_PI) * pow(t, z+0.5) * exp(-t) * x;
}


//The name for this macro can be at most 30 characters total
#define _class_print_species_(name,type) \
printf("-> %-30s Omega = %-15g , omega = %-15g\n",name,pba->Omega0_##type,pba->Omega0_##type*pba->h*pba->h);

/** list of possible types of spatial curvature */

enum spatial_curvature {flat,open,closed};

/** list of possible parametrisations of the DE equation of state */

enum equation_of_state {CLP,EDE};
enum ede_parametrization {tracker,pheno_axion};

/**
 * All background parameters and evolution that other modules need to know.
 *
 * Once initialized by the backgound_init(), contains all necessary
 * information on the background evolution (except thermodynamics),
 * and in particular, a table of all background quantities as a
 * function of time and scale factor, used for interpolation in other
 * modules.
 */
enum scf_pot{
  pol_times_exp, /** scf_potential set to pol_times_exp:V equals ((\phi-B)^\alpha + A)exp(-lambda*phi), see http://arxiv.org/abs/astro-ph/9908085.*/
  double_exp, /** scf_potential set to double_exp: V equals \Lambda_1^4e^{-\lambda\phi}+\Lambda_2^4e^{-\mu\phi} */
  axion, /** scf_potential set to axion: V equals m^2f^2(1-cos(phi/f)) */
  phi_2n, /** scf_potential set to axion: V equals V0((phi)^2n) */
  axionquad, /* scf_potential set to axion quadratic form: V = m^2phi^2/2 */
  ax_cos_cubed
};
struct background
{
  /** @name - input parameters initialized by user in input module
   *  (all other quantities are computed in this module, given these parameters
   *   and the content of the 'precision' structure)
   *
   * The background cosmological parameters listed here form a parameter
   * basis which is directly usable by the background module. Nothing
   * prevents from defining the input cosmological parameters
   * differently, and to pre-process them into this format, using the input
   * module (this might require iterative calls of background_init()
   * e.g. for dark energy or decaying dark matter). */

  //@{

  double H0; /**< \f$ H_0 \f$: Hubble parameter (in fact, [\f$H_0/c\f$]) in \f$ Mpc^{-1} \f$ */

  double Omega0_g; /**< \f$ \Omega_{0 \gamma} \f$: photons */

  double T_cmb; /**< \f$ T_{cmb} \f$: current CMB temperature in Kelvins */

  double Omega0_b; /**< \f$ \Omega_{0 b} \f$: baryons */

  double Omega0_cdm; /**< \f$ \Omega_{0 cdm} \f$: cold dark matter */

  double Omega0_lambda; /**< \f$ \Omega_{0_\Lambda} \f$: cosmological constant */

  double Omega0_fld; /**< \f$ \Omega_{0 de} \f$: fluid */

  enum equation_of_state fluid_equation_of_state; /**< parametrisation scheme for fluid equation of state */
  enum ede_parametrization ede_parametrization; /**< parametrisation scheme for ede type of fluid */

  double w0_fld; /**< \f$ w0_{DE} \f$: current fluid equation of state parameter */
  double wa_fld; /**< \f$ wa_{DE} \f$: fluid equation of state parameter derivative */
  double Omega_EDE; /**< \f$ wa_{DE} \f$: Early Dark Energy density parameter */

  double nu_fld; /**< parameter controlling width of ede transition from -1 to whatever w_f is */
  double n_pheno_axion; /**< exponent of (1-cos) in axion potential that fld is mimicing such that final eq of state of fld is n-1/n+1 */
  double Omega_fld_ac; /**< fractional energy density of EDE at a_c */
  double n_cap_infinity; /**< n_fld higher than this is assumed to be infinite. Helps set w_fld_final = 1 */
  double w_fld_f; /**< Final eq of state of EDE */
  double log10_a_c; /**< log10(a_c) to set a_c */
  double a_peak; /**< scale factor when EDE energy density peaks */
  double f_ede_peak; /**< maximum fractional energy density of EDE */
  short cs2_is_wn; /**< Switch in input to check if wn = cs2, read in cs2 and accordingly reset n_pheno_axion or w_final */
  double precision_newton_method_x; /**<precision at which we want the shooting to converge; default =1e-4 for theta_s */
  double precision_newton_method_F; /**<precision at which we want the shooting to converge; default =1e-6 for theta_s */
  double cs2_fld; /**< \f$ c^2_{s~DE} \f$: sound speed of the fluid
		     in the frame comoving with the fluid (so, this is
		     not [delta p/delta rho] in the synchronous or
		     newtonian gauge!) */

  short use_ppf; /**< flag switching on PPF perturbation equations
                    instead of true fluid equations for
                    perturbations. It could have been defined inside
                    perturbation structure, but we leave it here in
                    such way to have all fld parameters grouped. */

  double c_gamma_over_c_fld; /**< ppf parameter defined in eq. (16) of 0808.3125 [astro-ph] */

  double Omega0_ur; /**< \f$ \Omega_{0 \nu r} \f$: ultra-relativistic neutrinos */

  double Omega0_idr; /**< \f$ \Omega_{0 idr} \f$: interacting dark radiation */
  double T_idr;      /**< \f$ T_{idr} \f$: current temperature of interacting dark radiation in Kelvins */

  double Omega0_idm_dr; /**< \f$ \Omega_{0 idm_dr} \f$: dark matter interacting with dark radiation */

  double Omega0_dcdmdr; /**< \f$ \Omega_{0 dcdm}+\Omega_{0 dr} \f$: decaying cold dark matter (dcdm) decaying to dark radiation (dr) */

  double Gamma_dcdm; /**< \f$ \Gamma_{dcdm} \f$: decay constant for decaying cold dark matter */
  double tau_dcdm; /**< \f$ \tau_{dcdm} \f$: decay lifetime for decaying cold dark matter, in s^{-1} */


  double Omega_ini_dcdm;    /**< \f$ \Omega_{ini,dcdm} \f$: rescaled initial value for dcdm density (see 1407.2418 for definitions) */

  double Omega0_scf;        /**< \f$ \Omega_{0 scf} \f$: scalar field */
  short scf_evolve_as_fluid; /** set to false to only evolve KG equations, otherwise - switch to fluid when necessary. To be used in perturbation module*/
  double threshold_scf_fluid_m_over_H; /** if scf_evolve_as_fluid set to true, the scf will be modeled as a fluid once m/H drops below threshold_scf_fluid_m_over_H */
  double security_small_Omega_scf; /** enforce fluid when Om_scf is below  security_small_Omega_scf even if scf_evolve_as_fluid = False to avoid code crashing; harmless due to the smallness of Om_scf */
  short attractor_ic_scf;   /**< whether the scalar field has attractor initial conditions */
  double phi_ini_scf;       /**< \f$ \phi(t_0) \f$: scalar field initial value */
  double phi_prime_ini_scf; /**< \f$ d\phi(t_0)/d\tau \f$: scalar field initial derivative wrt conformal time */
  enum scf_pot scf_potential; /**< List of currently implement potential for a scalar field */
  double * scf_parameters;  /**< list of parameters describing the scalar field potential */
  int scf_parameters_size;  /**< size of scf_parameters */
  int scf_tuning_index;     /**< index in scf_parameters used for tuning */
  double m_scf;
  double f_axion;
  double alpha_squared;
  double power_of_mu;
  double log10_f_axion;
  double log10_m_axion;
  double log10_axion_ac;
  double V0_phi2n;
  double a_c;
  double log10_fraction_axion_ac;
  double adptative_stepsize;
  double omega_axion;
  double Omega0_axion;
  double Omega_axion_ac;
  double log10_z_c; //
  double zc_is_zeq; //enforce zc=zeq, with zeq=(Omega0_b+Omega0_cdm)/(Omega0_g+Omega0_ur);
  double f_ede; // TK added doubles to fill with values of the exact z_c and fraction_ede eventually
  double phi_scf_c; // Added for debugging. Trying to see whether the value of phi at z_c is really 7/8 phi_ini
  double n_axion;
  double w_scf;
  double cs2_scf;

  double n_axion_security;
  short scf_kg_eq;    /**< evolve scalar field with KG equations */
  short scf_fluid_eq;    /**< evolve scalar field with KG equations */
  short scf_evolve_like_axionCAMB; /**< evolve scalar field perturbations like axionCAMB */
  short scf_has_perturbations; /** do scalar field perts */
  short loop_over_background_for_closure_relation; /** do we want to loop over background?*/
  double precision_loop_over_background;
  //double scf_lambda; /**< \f$ \lambda \f$ : scalar field exponential potential slope */
  //double scf_alpha;  /**< \f$ \alpha \f$ : Albrecht-Skordis polynomial slope */
  //double scf_B; /**< \f$ \alpha \f$ : Albrecht-Skordis field shift */
  //double scf_A; /**< \f$ \alpha \f$ : Albrecht-Skordis offset */

  double Omega0_k; /**< \f$ \Omega_{0_k} \f$: curvature contribution */

  int N_ncdm;                            /**< Number of distinguishable ncdm species */
  double * M_ncdm;                       /**< vector of masses of non-cold relic:
                                             dimensionless ratios m_ncdm/T_ncdm */
  double * Omega0_ncdm, Omega0_ncdm_tot; /**< Omega0_ncdm for each species and for the total Omega0_ncdm */
  double * deg_ncdm, deg_ncdm_default;   /**< vector of degeneracy parameters in factor
                                             of p-s-d: 1 for one family of neutrinos
                                             (= one neutrino plus its anti-neutrino,
                                             total g*=1+1=2, so deg = 0.5 g*); and its
					     default value */

  /* the following parameters help to define the analytical ncdm phase space distributions (p-s-d) */
  double * T_ncdm,T_ncdm_default;       /**< list of 1st parameters in
					     p-s-d of non-cold relics:
					     relative temperature
					     T_ncdm1/T_gamma; and its
					     default value */
  double * ksi_ncdm, ksi_ncdm_default;  /**< list of 2nd parameters in
					     p-s-d of non-cold relics:
					     relative chemical potential
					     ksi_ncdm1/T_ncdm1; and its
					     default value */
  double * ncdm_psd_parameters;         /**< list of parameters for specifying/modifying
                                             ncdm p.s.d.'s, to be customized for given model
                                             (could be e.g. mixing angles) */
  /* end of parameters for analytical ncdm p-s-d */

  /* the following parameters help to define tabulated ncdm p-s-d passed in file */
  int * got_files;                      /**< list of flags for each species, set to true if
					     p-s-d is passed through file */
  char * ncdm_psd_files;                /**< list of filenames for tabulated p-s-d */
  /* end of parameters for tabulated ncdm p-s-d */

  //@}

  /** @name - related parameters */

  //@{

  double h; /**< reduced Hubble parameter */
  double age; /**< age in Gyears */
  double conformal_age; /**< conformal age in Mpc */
  double K; /**< \f$ K \f$: Curvature parameter \f$ K=-\Omega0_k*a_{today}^2*H_0^2\f$; */
  int sgnK; /**< K/|K|: -1, 0 or 1 */
  double * m_ncdm_in_eV; /**< list of ncdm masses in eV (inferred from M_ncdm and other parameters above) */
  double Neff; /**< so-called "effective neutrino number", computed at earliest time in interpolation table */
  double Omega0_dcdm; /**< \f$ \Omega_{0 dcdm} \f$: decaying cold dark matter */
  double Omega0_dr; /**< \f$ \Omega_{0 dr} \f$: decay radiation */
  double Omega0_m;  /**< total non-relativistic matter today */
  double Omega0_r;  /**< total ultra-relativistic radiation today */
  double Omega0_de; /**< total dark energy density today, currently defined as 1 - Omega0_m - Omega0_r - Omega0_k */
  double a_eq;      /**< scale factor at radiation/matter equality */
  double H_eq;      /**< Hubble rate at radiation/matter equality [Mpc^-1] */
  double z_eq;      /**< redshift at radiation/matter equality */
  double tau_eq;    /**< conformal time at radiation/matter equality [Mpc] */

  //@}

  /** @name - other background parameters */

  //@{

  double a_today; /**< scale factor today (arbitrary and irrelevant for most purposes) */

  //@}

  /** @name - all indices for the vector of background (=bg) quantities stored in table */

  //@{

  int index_bg_a;             /**< scale factor */
  int index_bg_H;             /**< Hubble parameter in \f$Mpc^{-1}\f$ */
  int index_bg_H_prime;       /**< its derivative w.r.t. conformal time */

  /* end of vector in short format, now quantities in normal format */

  int index_bg_rho_g;         /**< photon density */
  int index_bg_rho_b;         /**< baryon density */
  int index_bg_rho_cdm;       /**< cdm density */
  int index_bg_rho_lambda;    /**< cosmological constant density */
  int index_bg_rho_fld;       /**< fluid density */
  int index_bg_w_fld;         /**< fluid equation of state */
  int index_bg_Omega_fld;     /**< fluid fractional energy density Omega_fld / Omega_tot as a function of a, needed to calculate the "peak" of pheno_axion EDE type fluid */
  int index_bg_rho_ur;        /**< relativistic neutrinos/relics density */
  int index_bg_rho_idm_dr;    /**< density of dark matter interacting with dark radiation */
  int index_bg_rho_idr;       /**< density of interacting dark radiation */
  int index_bg_rho_dcdm;      /**< dcdm density */
  int index_bg_rho_dr;        /**< dr density */

  int index_bg_phi_scf;       /**< scalar field value */
  int index_bg_phi_prime_scf; /**< scalar field derivative wrt conformal time */
  int index_bg_V_scf;         /**< scalar field potential V */
  int index_bg_dV_scf;        /**< scalar field potential derivative V' */
  int index_bg_ddV_scf;       /**< scalar field potential second derivative V'' */
  int index_bg_rho_scf;       /**< scalar field energy density */
  int index_bg_Omega_scf;       /**< scalar field fractional energy density */
  int index_bg_p_scf;         /**< scalar field pressure */
  int index_bg_p_prime_scf;         /**< scalar field pressure */
  int index_bg_w_scf;         /**< scalar field e.o.s. */
  int index_bg_dw_scf;         /**< scalar field derivative of e.o.s. w/r to tau */
  int index_bg_ddw_scf;         /**< scalar field double derivative of e.o.s. w/r to tau*/
  int index_bg_rho_ncdm1;     /**< density of first ncdm species (others contiguous) */
  int index_bg_p_ncdm1;       /**< pressure of first ncdm species (others contiguous) */
  int index_bg_pseudo_p_ncdm1;/**< another statistical momentum useful in ncdma approximation */

  int index_bg_rho_tot;       /**< Total density */
  int index_bg_p_tot;         /**< Total pressure */
  int index_bg_p_tot_prime;   /**< Conf. time derivative of total pressure */

  int index_bg_Omega_r;       /**< relativistic density fraction (\f$ \Omega_{\gamma} + \Omega_{\nu r} \f$) */

  /* end of vector in normal format, now quantities in long format */

  int index_bg_rho_crit;      /**< critical density */
  int index_bg_Omega_m;       /**< non-relativistic density fraction (\f$ \Omega_b + \Omega_cdm + \Omega_{\nu nr} \f$) */
  int index_bg_conf_distance; /**< conformal distance (from us) in Mpc */
  int index_bg_ang_distance;  /**< angular diameter distance in Mpc */
  int index_bg_lum_distance;  /**< luminosity distance in Mpc */
  int index_bg_time;          /**< proper (cosmological) time in Mpc */
  int index_bg_rs;            /**< comoving sound horizon in Mpc */

  int index_bg_D;             /**< scale independent growth factor D(a) for CDM perturbations */
  int index_bg_f;             /**< corresponding velocity growth factor [dlnD]/[dln a] */

  int bg_size_short;  /**< size of background vector in the "short format" */
  int bg_size_normal; /**< size of background vector in the "normal format" */
  int bg_size;        /**< size of background vector in the "long format" */

  //@}

  /** @name - background interpolation tables */

  //@{

  int bt_size;               /**< number of lines (i.e. time-steps) in the array */
  double * tau_table;        /**< vector tau_table[index_tau] with values of \f$ \tau \f$ (conformal time) */
  double * z_table;          /**< vector z_table[index_tau] with values of \f$ z \f$ (redshift) */
  double * background_table; /**< table background_table[index_tau*pba->bg_size+pba->index_bg] with all other quantities (array of size bg_size*bt_size) **/

  //@}

  /** @name - table of their second derivatives, used for spline interpolation */

  //@{

  double * d2tau_dz2_table; /**< vector d2tau_dz2_table[index_tau] with values of \f$ d^2 \tau / dz^2 \f$ (conformal time) */
  double * d2background_dtau2_table; /**< table d2background_dtau2_table[index_tau*pba->bg_size+pba->index_bg] with values of \f$ d^2 b_i / d\tau^2 \f$ (conformal time) */

  //@}


  /** @name - all indices for the vector of background quantities to be integrated (=bi)
   *
   * Most background quantities can be immediately inferred from the
   * scale factor. Only few of them require an integration with
   * respect to conformal time (in the minimal case, only one quantity needs to
   * be integrated with time: the scale factor, using the Friedmann
   * equation). These indices refer to the vector of
   * quantities to be integrated with time.
   * {B} quantities are needed by background_functions() while {C} quantities are not.
   */

  //@{

  int index_bi_a;       /**< {B} scale factor */
  int index_bi_rho_dcdm;/**< {B} dcdm density */
  int index_bi_rho_dr;  /**< {B} dr density */
  int index_bi_rho_fld; /**< {B} fluid density */
  int index_bi_rho_scf; /**< {B} scf density */
  int index_bi_phi_scf;       /**< {B} scalar field value */
  int index_bi_phi_prime_scf; /**< {B} scalar field derivative wrt conformal time */

  int index_bi_time;    /**< {C} proper (cosmological) time in Mpc */
  int index_bi_rs;      /**< {C} sound horizon */
  int index_bi_tau;     /**< {C} conformal time in Mpc */
  int index_bi_D;       /**< {C} scale independent growth factor D(a) for CDM perturbations. */
  int index_bi_D_prime; /**< {C} D satisfies \f$ [D''(\tau)=-aHD'(\tau)+3/2 a^2 \rho_M D(\tau) \f$ */

  int bi_B_size;        /**< Number of {B} parameters */
  int bi_size;          /**< Number of {B}+{C} parameters */

  //@}

  /** @name - flags describing the absence or presence of cosmological
      ingredients
      *
      * having one of these flag set to zero allows to skip the
      * corresponding contributions, instead of adding null contributions.
      */


  //@{

  short has_cdm;       /**< presence of cold dark matter? */
  short has_dcdm;      /**< presence of decaying cold dark matter? */
  short has_dr;        /**< presence of relativistic decay radiation? */
  short has_scf;       /**< presence of a scalar field? */
  short has_ncdm;      /**< presence of non-cold dark matter? */
  short has_lambda;    /**< presence of cosmological constant? */
  short has_fld;       /**< presence of fluid with constant w and cs2? */
  short has_ur;        /**< presence of ultra-relativistic neutrinos/relics? */
  short has_idr;       /**< presence of interacting dark radiation? */
  short has_idm_dr;    /**< presence of dark matter interacting with dark radiation? */
  short has_curvature; /**< presence of global spatial curvature? */

  //@}

  /**
   *@name - arrays related to sampling and integration of ncdm phase space distributions
   */


  //@{
  int * ncdm_quadrature_strategy; /**< Vector of integers according to quadrature strategy. */
  int * ncdm_input_q_size; /**< Vector of numbers of q bins */
  double * ncdm_qmax;   /**< Vector of maximum value of q */
  double ** q_ncdm_bg;  /**< Pointers to vectors of background sampling in q */
  double ** w_ncdm_bg;  /**< Pointers to vectors of corresponding quadrature weights w */
  double ** q_ncdm;     /**< Pointers to vectors of perturbation sampling in q */
  double ** w_ncdm;     /**< Pointers to vectors of corresponding quadrature weights w */
  double ** dlnf0_dlnq_ncdm; /**< Pointers to vectors of logarithmic derivatives of p-s-d */
  int * q_size_ncdm_bg; /**< Size of the q_ncdm_bg arrays */
  int * q_size_ncdm;    /**< Size of the q_ncdm arrays */
  double * factor_ncdm; /**< List of normalization factors for calculating energy density etc.*/

  //@}

  /**
   *@name - some flags needed for calling background functions
   */

  //@{

  short short_info;  /**< flag for calling background_at_eta and return little information */
  short normal_info; /**< flag for calling background_at_eta and return medium information */
  short long_info;   /**< flag for calling background_at_eta and return all information */

  short inter_normal;  /**< flag for calling background_at_eta and find position in interpolation table normally */
  short inter_closeby; /**< flag for calling background_at_eta and find position in interpolation table starting from previous position in previous call */

  //@}

  /** @name - technical parameters */

  //@{

  short shooting_failed;  /**< flag is set to true if shooting failed. */

  ErrorMsg shooting_error; /**< Error message from shooting failed. */

  short background_verbose; /**< flag regulating the amount of information sent to standard output (none if set to zero) */

  ErrorMsg error_message; /**< zone for writing error messages */

  //@}
};

/**
 * temporary parameters and workspace passed to the background_derivs function
 */

struct background_parameters_and_workspace {

  /* structures containing fixed input parameters (indices, ...) */
  struct background * pba;

  /* workspace */
  double * pvecback;

};

/**
 * temporary parameters and workspace passed to phase space distribution function
 */

struct background_parameters_for_distributions {

  /* structures containing fixed input parameters (indices, ...) */
  struct background * pba;

  /* Additional parameters */

  /* Index of current distribution function */
  int n_ncdm;

  /* Used for interpolating in file of tabulated p-s-d: */
  int tablesize;
  double *q;
  double *f0;
  double *d2f0;
  int last_index;

};

/**************************************************************/
/* @cond INCLUDE_WITH_DOXYGEN */
/*
 * Boilerplate for C++
 */
#ifdef __cplusplus
extern "C" {
#endif

  int background_at_tau(
			struct background *pba,
			double tau,
			short return_format,
			short inter_mode,
			int * last_index,
			double * pvecback
			);

  int background_tau_of_z(
                          struct background *pba,
                          double z,
                          double * tau
                          );

  int background_functions(
			   struct background *pba,
			   double * pvecback_B,
			   short return_format,
			   double * pvecback
			   );

  int background_w_fld(
                       struct background * pba,
                       double a,
                       double * w_fld,
                       double * dw_over_da_fld,
                       double * integral_fld);

  int background_init(
		      struct precision *ppr,
		      struct background *pba
		      );

  int background_free(
		      struct background *pba
		      );

  int background_free_input(
                            struct background *pba
                            );

  int background_free_noinput(
                    struct background *pba
                    );

  int background_indices(
			 struct background *pba
			 );

  int background_ncdm_distribution(
				  void *pba,
				  double q,
				  double * f0
				  );

  int background_ncdm_test_function(
				     void *pba,
				     double q,
				     double * test
				     );

  int background_ncdm_init(
			    struct precision *ppr,
			    struct background *pba
			    );


  int background_ncdm_momenta(
                             double * qvec,
                             double * wvec,
                             int qsize,
                             double M,
                             double factor,
                             double z,
                             double * n,
		             double * rho,
                             double * p,
                             double * drho_dM,
			     double * pseudo_p
                             );

  int background_ncdm_M_from_Omega(
				    struct precision *ppr,
				    struct background *pba,
					int species
				    );

  int background_solve(
		       struct precision *ppr,
		       struct background *pba
		       );

  int background_initial_conditions(
				    struct precision *ppr,
				    struct background *pba,
				    double * pvecback,
				    double * pvecback_integration
				    );

  int background_find_equality(
                               struct precision *ppr,
                               struct background *pba
                               );

  int background_output_titles(struct background * pba,
                               char titles[_MAXTITLESTRINGLENGTH_]
                               );

  int background_output_data(
                           struct background *pba,
                           int number_of_titles,
                           double *data);

  int background_derivs(
			 double z,
			 double * y,
			 double * dy,
			 void * parameters_and_workspace,
			 ErrorMsg error_message
			 );

  /** Scalar field potential and its derivatives **/
  double V_scf(
               struct background *pba,
               double phi
               );

  double dV_scf(
		struct background *pba,
		double phi
		);

  double ddV_scf(
                 struct background *pba,
                 double phi
                 );

  /** Coupling between scalar field and matter **/
  double Q_scf(
               struct background *pba,
               double phi,
               double phi_prime
               );

  /** Budget equation output */
  int background_output_budget(
               struct background* pba
               );

#ifdef __cplusplus
}
#endif

/**************************************************************/

/**
 * @name Some conversion factors and fundamental constants needed by background module:
 */

//@{

#define _Mpc_over_m_ 3.085677581282e22  /**< conversion factor from meters to megaparsecs */
/* remark: CAMB uses 3.085678e22: good to know if you want to compare  with high accuracy */

#define _Gyr_over_Mpc_ 3.06601394e2 /**< conversion factor from megaparsecs to gigayears
				         (c=1 units, Julian years of 365.25 days) */
#define _c_ 2.99792458e8            /**< c in m/s */
#define _G_ 6.67428e-11             /**< Newton constant in m^3/Kg/s^2 */
#define _eV_ 1.602176487e-19        /**< 1 eV expressed in J */
#define _eV_over_Kelvin_ 8.626e-5   /**< kB in eV/K */
#define _eV_over_joules_ 6.24150647996e+18 /**< eV/J */
#define _eV_over_ergs_ 6.24150647996e11 /**< eV/J */
#define _joules_over_ergs_ 1e-7 /**< eV/J */
#define _Sun_mass_over_kg_ 1.98855e30 /**< M_sun in kg */

/* parameters entering in Stefan-Boltzmann constant sigma_B */
#define _k_B_ 1.3806504e-23
#define _h_P_ 6.62606896e-34
/* remark: sigma_B = 2 pi^5 k_B^4 / (15h^3c^2) = 5.670400e-8
                   = Stefan-Boltzmann constant in W/m^2/K^4 = Kg/K^4/s^3 */

//@}

/**
 * @name Some limits on possible background parameters
 */

//@{

#define _H0_BIG_ 2./2997.9     /**< maximal \f$ H_0 \f$ in \f$ Mpc^{-1} (h=1.0) \f$ */
#define _H0_SMALL_ 0.3/2997.9  /**< minimal \f$ H_0 \f$ in \f$ Mpc^{-1} (h=0.3) \f$ */
#define _TCMB_BIG_ 3.5         /**< maximal \f$ T_{cmb} \f$ in K */
#define _TCMB_SMALL_ 2.       /**< minimal \f$ T_{cmb}  \f$ in K */
#define _TOLERANCE_ON_CURVATURE_ 1.e-5 /**< if \f$ | \Omega_k | \f$ smaller than this, considered as flat */
#define _OMEGAK_BIG_ 0.5     /**< maximal \f$ Omega_k \f$ */
#define _OMEGAK_SMALL_ -0.5  /**< minimal \f$ Omega_k \f$ */
#define _h_BIG_ 1.5            /**< maximal \f$ h \f$ */
#define _h_SMALL_ 0.3         /**< minimal \f$ h \f$ */
#define _omegab_BIG_ 0.039    /**< maximal \f$ omega_b \f$ */
#define _omegab_SMALL_ 0.005  /**< minimal \f$ omega_b \f$ */

//@}

/**
 * @name Some limits imposed in other parts of the module:
 */

//@{

#define _SCALE_BACK_ 0.1  /**< logarithmic step used when searching
			     for an initial scale factor at which ncdm
			     are still relativistic */

#define _PSD_DERIVATIVE_EXP_MIN_ -30 /**< for ncdm, for accurate computation of dlnf0/dlnq, q step is varied in range specified by these parameters */
#define _PSD_DERIVATIVE_EXP_MAX_ 2  /**< for ncdm, for accurate computation of dlnf0/dlnq, q step is varied in range specified by these parameters */

#define _zeta3_ 1.2020569031595942853997381615114499907649862923404988817922 /**< for quandrature test function */
#define _zeta5_ 1.0369277551433699263313654864570341680570809195019128119741 /**< for quandrature test function */

//@}


#endif
/* @endcond */

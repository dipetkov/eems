
#include "eems.hpp"

#include <boost/config.hpp>
#include <boost/program_options.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/detail/config_file.hpp>
namespace po = boost::program_options;

// The distance metric is a global variable, so that the pairwise_distance function can see it
// Choose 'euclidean' (the default) or 'greatcirc' (great circle distance)
string dist_metric;

int main(int argc, char** argv)
{
  // random seed
  long seed = time(NULL);
  string config_file;
  Params params;

  try {
    
    // Declare options allowed only on the command line
    po::options_description generic_options("Generic options");
    generic_options.add_options()
      ("help", "Produce this help message")
      ("params", po::value<string>(&config_file)->required(),
       "Specify the EEMS configuration file")
      ;
    
    // Declare options allowed on both the command line and in a config file
    po::options_description eems_options("EEMS options");
    eems_options.add_options()
      ("seed", po::value<long>(&seed),
       "Set the random seed")
      ("datapath", po::value<string>(&params.datapath),
       "Full path to a set of three files: datapath.coord, datapath.sites and datapath.outer.")
      ("mcmcpath", po::value<string>(&params.mcmcpath),
       "Full path to an output directory with write permission. Will be created if necessary.")
      ("prevpath", po::value<string>(&params.prevpath)->default_value(""),
       "Full path to previous output directory, i.e., the mcmcpath in a previous EEMS run.")
      ("gridpath", po::value<string>(&params.gridpath)->default_value(""),
       "Full path to a set of three files: gridpath.demes, gridpath.edges and gridpath.ipmap.")
      // Interpret the input nDemes as a suggestion because there are some further constraints
      // when constructing a regular triangular grid in an irregular habitat
      ("nDemes", po::value<int>(&params.nDemes),
       "Number of demes, roughly. EEMS constructs a regular triangular grid with ~nDemes vertices.")
      ("nIndiv", po::value<int>(&params.nIndiv),
       "Number of samples. Should match the number of rows in datapath.sites.")
      ("nSites", po::value<int>(&params.nSites),
       "Number of microsatellites. Should match the number of columns in datapath.sites.")
      ("diploid", po::value<bool>(&params.diploid)->default_value(true),
       "Logical. Indicates whether the species is diploid (true) or haploid (false).")
      ("distance", po::value<string>(&dist_metric)->default_value("euclidean"),
       "Distance metric. Either 'euclidean' or 'greatcirc', i.e., great circle.")
      ("numMCMCIter", po::value<int>(&params.numMCMCIter),
       "Number of Markov Chain Monte Carlo iterations.")
      ("numBurnIter", po::value<int>(&params.numBurnIter),
       "Number of burn-in iterations to be discarded before sampling from posterior.")
      ("numThinIter", po::value<int>(&params.numThinIter),
       "Number of thinning iterations to be discarded between sampling from posterior.")
      ("mSeedsProposalS2", po::value<double>(&params.mSeedsProposalS2)->default_value(0.01),
       "Variance of normal proposals to update the seeds of the migration tiles.")
      ("qSeedsProposalS2", po::value<double>(&params.qSeedsProposalS2)->default_value(0.1),
       "Variance of normal proposals to update the seeds of the diversity tiles.")
      ("mEffctProposalS2", po::value<double>(&params.mEffctProposalS2)->default_value(0.1),
       "Variance of normal proposals to update the log10 rates of the migration tiles.")
      ("qEffctProposalS2", po::value<double>(&params.qEffctProposalS2)->default_value(0.001),
       "Variance of normal proposals to update the log10 rates of the diversity tiles.")
      ("mrateMuProposalS2", po::value<double>(&params.mrateMuProposalS2)->default_value(0.01),
       "Variance of normal proposals to update the overall mean migration rate, on the log10 scale.")
      ("qVoronoiPr", po::value<double>(&params.qVoronoiPr)->default_value(0.25),
       "With prob. qVoronoiPr, update diversity Voronoi; with prob. 1-qVoronoiPr, update migration Voronoi.")
      ("qrateShape", po::value<double>(&params.qrateShape_2)->default_value(0.001),
       "Shape hyperparameter for the diversity rates variance, qrateS2 ~ invgamma(qrateShape, qrateScale)")
      ("mrateShape", po::value<double>(&params.mrateShape_2)->default_value(0.001),
       "Shape hyperparameter for the migration rates variance, mrateS2 ~ invgamma(mrateShape, mrateScale)")
      ("sigmaShape", po::value<double>(&params.sigmaShape_2)->default_value(0.001),
       "Shape hyperparameter for the scaling factor sigma^2 ~ invgamma(sigmaShape, sigmaScale)")
      ("qrateScale", po::value<double>(&params.qrateScale_2)->default_value(1.0),
       "Scale hyperparameter for the diversity rates variance, qrateS2 ~ invgamma(qrateShape, qrateScale)")
      ("mrateScale", po::value<double>(&params.mrateScale_2)->default_value(1.0),
       "Scale hyperparameter for the migration rates variance, mrateS2 ~ invgamma(mrateShape, mrateScale)")
      ("sigmaScale", po::value<double>(&params.sigmaScale_2)->default_value(1.0),
       "Scale hyperparameter for the scaling factor sigma^2 ~ invgamma(sigmaShape, sigmaScale)")
      ("negBiProb", po::value<double>(&params.negBiProb)->default_value(0.67),
       "Success probability for the number of Voronoi tiles ~ negbinom(negBiSize, negBiProb)")
      ("negBiSize", po::value<int>(&params.negBiSize)->default_value(10),
       "Size for the number of Voronoi tiles ~ negbinom(negBiSize, negBiProb)")
      ;
    
    po::options_description cmdline_options;
    cmdline_options.add(generic_options).add(eems_options);
    
    po::options_description config_file_options;
    config_file_options.add(eems_options);

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, cmdline_options), vm);
    po::notify(vm);
    
    if (vm.count("help")) {
      cout << cmdline_options << "\n";
      return 0;
    }
    
    ifstream instrm(config_file.c_str());
    if (instrm) {
      po::store(po::parse_config_file(instrm, config_file_options, true), vm);
      po::notify(vm);
    }

    params.seed = seed;
    params.check_input_arguments();
    
    EEMS eems(params);
    MCMC mcmc(params);
    
    boost::filesystem::path dir(eems.prevpath().c_str());
    if (exists(dir)) {
      cout << "Load final EEMS state from " << eems.prevpath() << endl << endl;
      eems.load_final_state();
    } else {
      cout << "Initialize EEMS random state" << endl << endl;
      eems.initialize_state();
    }
  
    if (!eems.start_eems(mcmc)) {
      cerr << "[RunEEMS] Error starting EEMS." << endl;
      return(EXIT_FAILURE);
    }
  
    Proposal proposal;

    while (!mcmc.finished) {
      
      switch ( eems.choose_move_type( ) ) {
      case Q_VORONOI_BIRTH_DEATH:
	eems.propose_birthdeath_qVoronoi(proposal);
	break;
      case M_VORONOI_BIRTH_DEATH:
	eems.propose_birthdeath_mVoronoi(proposal);
	break;
      case Q_VORONOI_POINT_MOVE:
	eems.propose_move_one_qtile(proposal);
	break;
      case M_VORONOI_POINT_MOVE:
	eems.propose_move_one_mtile(proposal);
	break;
      case Q_VORONOI_RATE_UPDATE:
	eems.propose_rate_one_qtile(proposal);
	break;
      case M_VORONOI_RATE_UPDATE:
	eems.propose_rate_one_mtile(proposal);
	break;
      case M_MEAN_RATE_UPDATE:
	eems.propose_overall_mrate(proposal);
	break;
      default:
	cerr << "[RunEEMS] Unknown move type" << endl;
	return(EXIT_FAILURE);
      }
	
      mcmc.add_to_total_moves(proposal.move);
      if (eems.accept_proposal(proposal)) { mcmc.add_to_okay_moves(proposal.move); }
      if (params.testing) { eems.check_ll_computation( ); }
      
      eems.update_sigma2( );
      eems.update_hyperparams( );
      mcmc.end_iteration( );
      
      // Check whether to save the current parameter state,
      // as the thinned out iterations are not saved
      if (mcmc.to_save_iteration()) {
	eems.print_iteration(mcmc);
	eems.save_iteration(mcmc);
      }      
    }
    
    eems.output_results(mcmc);
        
  } catch(exception& e) {
    cerr << e.what() << endl;
    return(EXIT_FAILURE);
  }    
  
  return(0);
}

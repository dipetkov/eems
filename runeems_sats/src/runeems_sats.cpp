
#include "eems.hpp"

// The distance metric is a global variable, so that
// the pairwise_distance function can see it
// Choose 'euclidean' (the default) or 'greatcirc' (great circle distance)
string dist_metric;

int main(int argc, char** argv)
{
  try {
    
    long seed_from_command_line = 1L;
    string params_file; bool error;
    
    po::options_description options("EEMS options from command line");
    po::variables_map vm;
    options.add_options()
      ("help", "Produce this help message")
      ("seed", po::value<long>(&seed_from_command_line)->default_value(time(NULL)), "Set the random seed")
      ("params", po::value<string>(&params_file)->required(), "Specify input parameter file") ;

    po::store(po::parse_command_line(argc, argv, options), vm);
    po::notify(vm);
  
    if(vm.count("help")) {
      cerr << options << endl; return(EXIT_FAILURE);
    }
    
    Params params(params_file,seed_from_command_line);
    error = params.check_input_params( );
    if (error) {
      cerr << "[RunEEMS] Error parametrizing EEMS." << endl;
      return(EXIT_FAILURE);      
    }
    
    // Specify the distance metric in the params.ini file
    dist_metric = params.distance;

    EEMS eems(params);
    MCMC mcmc(params);
    
    boost::filesystem::path dir(eems.prevpath().c_str());
    if (exists(dir)) {
      cerr << "Load final EEMS state from " << eems.prevpath() << endl << endl;
      eems.load_final_state();
    } else {
      cerr << "Initialize EEMS random state" << endl << endl;
      eems.initialize_state();
    }
    
    error = eems.start_eems(mcmc);
    if (error) {
      cerr << "[RunEEMS] Error starting EEMS." << endl;
      return(EXIT_FAILURE);
    }
    
    Proposal proposal;
    
    while (!mcmc.finished) {

      switch ( eems.choose_move_type( ) ) {
      case Q_VORONOI_BIRTH_DEATH:
	eems.birthdeath_qVoronoi(proposal);
	break;
      case M_VORONOI_BIRTH_DEATH:
	eems.birthdeath_mVoronoi(proposal);
	break;
      case Q_VORONOI_POINT_MOVE:
	eems.move_qVoronoi(proposal);
	break;
      case M_VORONOI_POINT_MOVE:
	eems.move_mVoronoi(proposal);
	break;
      case Q_VORONOI_RATE_UPDATE:
	eems.propose_qEffcts(proposal);
	break;
      case M_VORONOI_RATE_UPDATE:
	eems.propose_mEffcts(proposal);
	break;
      case M_MEAN_RATE_UPDATE:
	eems.propose_mrateMu(proposal);
	break;
      default:
	cerr << "[RunEEMS] Unknown move type" << endl;
	return(EXIT_FAILURE);
      }

      mcmc.add_to_total_moves(proposal.type);
      if (eems.accept_proposal(proposal)) { mcmc.add_to_okay_moves(proposal.type); }
      if (params.testing) { eems.check_ll_computation( ); }
      
      eems.update_sigma2( );
      eems.update_hyperparams( );
      mcmc.end_iteration( );
      
      // Check whether to save the current parameter state,
      // as the thinned out iterations are not saved
      int iter = mcmc.to_save_iteration( );
      if (iter>=0) {
	eems.print_iteration(mcmc);
	eems.save_iteration(mcmc);
      }      
    }
    
    error = eems.output_results(mcmc);
    if (error) { cerr << "[RunEEMS] Error saving eems results to " << eems.mcmcpath() << endl; }
    
  } catch(exception& e) {
    cerr << e.what() << endl;
    return(EXIT_FAILURE);
  }    
  
  return(0);
}

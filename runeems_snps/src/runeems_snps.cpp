
#include "util.hpp"
#include "eems.hpp"
#include "mcmc.hpp"

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
    
    cerr << "Input parameters: " << endl << params << endl
	 << "Initial log prior: " << eems.prior( ) << endl
	 << "Initial log llike: " << eems.likelihood( ) << endl << endl;
    
    Proposal proposal;
    
    while (!mcmc.finished) {
      
      if (!(mcmc.currIter%1000)) {
	cerr << "Iteration " << mcmc.currIter << "..." << endl;
      }
      
      mcmc.start_iteration( );
      eems.update_sigma2( );

      while (!mcmc.iterDone) {

	// There are 5 'steps' as EEMS cycles through 5 types of proposals:
	// birth/death, move a tile, update the rate of  a tile, update the
	// mean migration rate, update the degrees of freedom
	double u = eems.runif( );
	switch (mcmc.currStep) {
	case 0:
	  // Make a birth or death proposal, with equal probability
	  // Choose to make a birth/death proposal in the Voronoi
	  // tessellation of the diversity rates with probability
	  // params.qVoronoiPr (which is 0.05 by default). Otherwise,
	  // make a birth/death proposal in the Voronoi tessellation
	  // of the migration rates
	  if (u < params.qVoronoiPr) {
	    eems.birthdeath_qVoronoi(proposal);
	  } else {
	    eems.birthdeath_mVoronoi(proposal);
	  }
	  break;
	case 1:
	  // Propose to move an existing tile within the habitat
	  if (u < params.qVoronoiPr) {
	    eems.move_qVoronoi(proposal);
	  } else {
	    eems.move_mVoronoi(proposal);
	  }
	  break;
	case 2:
	  // Propose to update the rate parameter of an existing tile
	  if (u < params.qVoronoiPr) {
	    eems.propose_qEffcts(proposal);
	  } else {
	    eems.propose_mEffcts(proposal);
	  }
	  break;
	case 3:
	  // Propose to update the overall (log10) migration rate
	  eems.propose_mrateMu(proposal);
	  break;
	default:
	  // Propose to update the degrees of freedom
	  eems.propose_df(proposal,mcmc);
	}
	
	mcmc.add_to_total_moves(proposal.type);
	if (eems.accept_proposal(proposal)) { mcmc.add_to_okay_moves(proposal.type); }
	if (params.testing) { eems.check_ll_computation( ); }
	mcmc.change_update( );
      }
      
      eems.update_hyperparams( );
      
      // Check whether to save the current parameter state,
      // as the thinned out iterations are not saved
      int iter = mcmc.to_save_iteration( );
      if (iter>=0) {
	eems.print_iteration(mcmc);
	eems.save_iteration(mcmc);
      }      
      mcmc.end_iteration( );
    }
    
    error = eems.output_results(mcmc);
    if (error) { cerr << "[RunMCMC] Error saving eems results to " << eems.mcmcpath() << endl; }
    
    cerr << "Final log prior: " << eems.prior( ) << endl
	 << "Final log llike: " << eems.likelihood( ) << endl;
    
  } catch(exception& e) {
    cerr << e.what() << endl;
    return(EXIT_FAILURE);
  }    
  
  return(0);
}


#include "util.hpp"
#include "eems.hpp"
#include "mcmc.hpp"


int main(int argc, char** argv)
{

  long seed = 1L;
  string params_file;

  po::options_description options("Options");
  po::variables_map vm;
  options.add_options()
    ("help", "Produce this help message")
    ("seed",po::value<long>(&seed)->default_value(time(NULL)), "Set the random seed")
    ("params",po::value<string>(),"Specify input parameter file") ;

  po::store(po::parse_command_line(argc, argv, options), vm);
  po::notify(vm);
  
  if(vm.count("help")) {
    cerr << options << endl;
    return EXIT_FAILURE;
  }
  if(vm.count("params")) {
    params_file = vm["params"].as<string>();
  } else {
    cerr << "[EEMS::Params] Please provide a params file with the following information:" << endl
	 << "               datapath, mcmcpath, nIndiv, nSites, nDemes" << endl;
    return EXIT_FAILURE;
  }
  
  cerr << "[EEMS::Params] Random seed = " << seed << endl;
  
  Params params(params_file);
  EEMS eems(params,seed);
  MCMC mcmc(params);

  eems.initialize(mcmc);
  cerr << fixed << setprecision(2)
       << "[RunMCMC] Initial log prior = " << eems.eval_prior( ) << endl
       << "          Initial log llike = " << eems.eval_likelihood( ) << endl;

  Proposal proposal;

  while (!mcmc.finished) {

    if (!mod(mcmc.currIter,100)) {
      cerr << "Iteration " << mcmc.currIter << "..." << endl;
    }
    
    mcmc.start_iteration( );
    eems.update_s2loc( );

    while (!mcmc.iterDone) {

      switch (mcmc.currType) {
      case 0:
	eems.propose_qEffcts(proposal,mcmc);
	break;
      case 1:
	eems.move_qVoronoi(proposal,mcmc);
	break;
      case 2:
	eems.birthdeath_qVoronoi(proposal,mcmc);
	break;
      case 3:
	eems.propose_mEffcts(proposal,mcmc);
	break;
      case 4:
	eems.propose_mrateMu(proposal);
	break;
      case 5:
	eems.move_mVoronoi(proposal,mcmc);
	break;
      default:
	eems.birthdeath_mVoronoi(proposal,mcmc);
      }

      mcmc.addToTotalMoves( );
      if (eems.accept_proposal(proposal)) { mcmc.addToOkayMoves( ); }
      if (params.testing)           { eems.check_ll_computation( ); }
      mcmc.change_update(eems.num_qtiles(),eems.num_mtiles());
    }

    eems.update_hyperparams( );
    
    // Check whether to save the current parameter state
    int iter = mcmc.to_save_iteration( );
    if (iter>=0) {
      cerr << "Ending iteration " << mcmc.currIter << " with acceptance proportions:" << endl;    
      mcmc.output_proportions(cerr);
      eems.report_iteration(mcmc.currIter);
      eems.save_iteration(iter);
    }

    mcmc.end_iteration( );
  }

  cerr << fixed << setprecision(2)
       << "[RunMCMC] Final log prior = " << eems.eval_prior( ) << endl
       << "          Final log llike = " << eems.eval_likelihood( ) << endl;

  eems.output_results(mcmc);

  return 0;
}

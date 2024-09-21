#include "choice.hpp"

/*************************************************
*		Initial values							 *
**************************************************/
char*		Choice::inFileName = nullptr;
char*		Choice::outPath = nullptr;


Method_t 		Choice::method(Method_t::OPT);

double			Choice::tolerance(0.2);
unsigned int	Choice::budget(1E06);
unsigned int	Choice::trial(3);


const char short_choices[] = "i:l:s:t:u:";

static struct choice long_choices[] = {
	{"method", required_argument, 0, METHOD},
	{"budget", required_argument, 0, SIZE_BUDGET},
	{"tolerance", required_argument, 0, SIZE_TOLERANCE}, 
	{"trial", required_argument, 0, SIZE_TRIAL}, 
	{0, 0, 0, 0}
};

bool Choice::parse(int argc, char **argv)
{

	int opt = 0;
	int choice_index = 0;
	while(1)
	{
		opt = getopt_long(argc, argv, short_choices, long_choices, &choice_index);
		if (-1 == opt)
			break;

		switch(opt)
		{
			case SIZE_BUDGET:
				Choice::budget = atoi(optarg);
				break;

			case SIZE_TOLERANCE:
				Choice::tolerance = atof(optarg);
				break;

			case SIZE_TRIAL:
				Choice::trial = atoi(optarg);
				break;

			case METHOD:
				if (strcasecmp(optarg, "naive") == 0)
					Choice::method = Method_t::NAIVE;
				else if (strcasecmp(optarg, "simple") == 0)
					Choice::method = Method_t::HASH;
				else
				{
					cout << "Unrecognized method. OPT method is used by default" << endl;
					Choice::method = Method_t::OPT;
				}
				break;

			default:
				if (long_choices[choice_index].flag != 0)
					break;

				cerr << "Invalid parameter " << argv[optind] << " " << optarg << endl;
				exit(1);
		}
	}

	if (optind < argc - 1)
	{
		Choice::inFileName = argv[optind++];
		Choice::outPath = argv[optind++];
	}
	else
	{
		cout << "Invalid parameter" << endl;
		exit(1);
	}

	return true;

}

void Choice::print()
{

	cout << Choice::inFileName << "\t" << Choice::outPath << "\t" << Choice::budget << "\t" << Choice::tolerance << "\t" << Choice::trial << endl;

	cout << "Method: ";
	switch (Choice::method)
	{
		case Method_t::NAIVE:
			cout << "Naive";
			break;
		case Method_t::HASH:
			cout << "Simple";
			break;
		case Method_t::OPT:
			cout << "OPT";
			break;
	}
	cout << endl;


	return;
}


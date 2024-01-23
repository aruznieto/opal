/***************************************************************/
//
//Copyright (c) 2019 Esteban Egea-Lopez http://girtel.upct.es/~simulations
//
/**************************************************************/

#include "timer.h"
#include "Opal.h"
#include "configuration.h"
#include <memory>
using namespace opal;
using namespace optix;
//Tests. Compile as exe
#ifdef _WIN32
//Define wgetopt as a substitute for getopt in Windows

int     opterr = 1,             /* if error message should be printed */
optind = 1,             /* index into parent argv vector */
optopt,                 /* character checked for validity */
optreset;               /* reset getopt */
char* optarg;                /* argument associated with option */

#define BADCH   (int)'?'
#define BADARG  (int)':'
#define EMSG    ""

/*
 * getopt --
 *      Parse argc/argv argument vector.
 */
int wgetopt(int nargc, char* const nargv[], const char* ostr)
{
	static char* place = EMSG;              /* option letter processing */
	const char* oli;                              /* option letter list index */

	if (optreset || !*place) {              /* update scanning pointer */
		optreset = 0;
		if (optind >= nargc || *(place = nargv[optind]) != '-') {
			place = EMSG;
			return (-1);
		}
		if (place[1] && *++place == '-') {      /* found "--" */
			++optind;
			place = EMSG;
			return (-1);
		}
	}                                       /* option letter okay? */
	if ((optopt = (int)*place++) == (int)':' ||
		!(oli = strchr(ostr, optopt))) {
		/*
		* if the user didn't specify '-' as an option,
		* assume it means -1.
		*/
		if (optopt == (int)'-')
			return (-1);
		if (!*place)
			++optind;
		if (opterr && *ostr != ':')
			(void)printf("illegal option -- %c\n", optopt);
		return (BADCH);
	}
	if (*++oli != ':') {                    /* don't need argument */
		optarg = NULL;
		if (!*place)
			++optind;
	}
	else {                                  /* need an argument */
		if (*place)                     /* no white space */
			optarg = place;
		else if (nargc <= ++optind) {   /* no arg */
			place = EMSG;
			if (*ostr == ':')
				return (BADARG);
			if (opterr)
				(void)printf("option requires an argument -- %c\n", optopt);
			return (BADCH);
		}
		else                            /* white space */
			optarg = nargv[optind];
		place = EMSG;
		++optind;
	}
	return (optopt);                        /* dump back option letter */
}

#else 
#include <unistd.h>
#endif

//Include the file with the tests you want to perform
#include "tests/depolarization.h"
#include "tests/tunnels.h"
#include "tests/tests.h"
#include "tests/rdn.h"
#include "tests/curvature.h"
#include "tests/dudley.h"
#include "tests/anvers.h"
#include "tests/anversDoa.h"
#include "tests/anversMimo.h"
#include "tests/rouxAP.h"
#include "tests/diffraction.h"
#include "tests/antennaGain.h"
#include "tests/cartagenaGSM.h"
#include "tests/parking.h"
#include "tests/lora.h"
#include "tests/lille.h"
#include "tests/exposure.h"
int main(int argc, char** argv)
{
	try {

		//Init log
		//Remove some info from the log
		loguru::g_preamble_date=false;
		loguru::g_preamble_time=false;
		loguru::g_preamble_thread=false;
		loguru::init(argc,argv);
		
		std::cout << "Running Opal with:  " ;
		for(int i = 0; i < argc; ++i) {
 		       std::cout << argv[i] << ' ';
		}
		std::cout<<std::endl;
		//float freq = 900e6f;
		//float freq = 5.9e9f;
		//float freq = 1.0e9f;

		//Some default values
		uint maxReflections=10u;
		bool printEnabled=false;
		bool subSteps=false;
		bool useExactSpeedOfLight=true;
		bool useDepolarization=true;
		bool usePenetration = false;
		bool useMultiGPU=true;
		bool useFastMath=true;
		float radius=0.04f;
        	ConfigurationParser* configuration = nullptr;

		std::string usage="Usage: opal [-options] \n  -r Max reflections E \n -p Enable OptiX rtPrintf on device to debug \n -s Use decimal degrees in angular spacing \n -c Use c=3e8 m/s. Default is c=299 792 458 m/s\n -d Enable depolarization \n -a Enable penetration \n -g Disable multiGPU \n -m disable fast_math \n -h Show help \n -t string to select a particular test (dependent on the program)\n -o Path to output folder\n -u Radius of the receiver sphere\n -f Path to JSON file with configuration keys\n";

		int c;
		int nr;
		std::string test;
		std::string outputPath;
#ifdef _WIN32
		while ((c = wgetopt (argc, argv, "r:pscdagmu:ht:o:")) != -1) {
#else 

		while ((c = getopt (argc, argv, "r:pscdagmu:ht:o:f:")) != -1) {
#endif
			switch (c) {
				case 'r':
					//std::cout<<optarg<<std::endl;
					nr=atoi(optarg);
					if (nr>=0) {

						maxReflections=nr;
					} else {
						fprintf(stderr, usage.c_str(),
								argv[0]);
						exit(EXIT_FAILURE);
					}
					break;
				case 'f':
                    			configuration = new ConfigurationParser(optarg);
					break;
				case 'u':
					//std::cout<<"-u="<<optarg<<std::endl;
					radius=atof(optarg);
					break;
				case 't':
					//std::cout<<optarg<<std::endl;
					test=optarg;
					break;
				case 'o':
					//std::cout<<optarg<<std::endl;
					outputPath=optarg;
					break;
				case 'p':
					printEnabled=true;
					break;
				case 's':
					subSteps=true;
					break;
				case 'c':
					useExactSpeedOfLight=false;
					break;
				case 'd':
					useDepolarization=true;
					break;
				case 'a':
					usePenetration=true;
					break;
				case 'g':
					useMultiGPU=false;
					break;
				case 'm':
					useFastMath=false;
					break;
				case 'h':
					std::cout<<usage<<std::endl;
					exit(0);
					break;

				default: /* '?' */
					fprintf(stderr, usage.c_str(),
							argv[0]);
					exit(EXIT_FAILURE);


			}
		}
//#endif
		

		std::cout<<"verbosity"<<loguru::g_stderr_verbosity<<std::endl;	
		OpalSceneManager* sceneManager=new OpalSceneManager();


		//Now set desired common features
		//Use of a configuration file probably will override this settings in the particular test
		if (!useExactSpeedOfLight) {
			sceneManager->useApproximateSpeedLight();
		}

		if (usePenetration) {
			std::cout<<"WARNING: penetration requires further test and works only with flat mesh and depolarization (no RDN or curved surfaces)" <<std::endl;		
			sceneManager->enablePenetration();
		}
		
		if (useDepolarization) {	
			sceneManager->enableDepolarization();
		}
		if (!useMultiGPU) {
			sceneManager->disableMultiGPU();	
		}
		if (!useFastMath) {
			sceneManager->disableFastMath();
		}
		sceneManager->setMaxReflections(maxReflections);

		//Set this for detailed debug of scene		
		//sceneManager->setUsageReport();

		//Enable debug of optix rays
		if (printEnabled) {
			sceneManager->setPrintEnabled(1024* 1024 * 1024);	
		}

		//Additional features are set in the particular tests and initContext is called then 
		

		
		

		
		//Finally, run the test


//Basic tests
		BasicTests basic(sceneManager,radius,useDepolarization);
		if (configuration) {
			basic.setConfiguration(configuration);
		}
		//basic.planeTest(2);
		//basic.useCallbacks(2);
		//basic.planeTestMultichannel(0);
	//	basic.freeSpace();
		//basic.quadTest(false,false);
		//basic.addCompoundDynamicMeshes();

		//basic.loadScenario();
		//basic.crossingTest(false,subSteps,useDepolarization,radius);
		//basic.crossingTestEfficient(false,subSteps,useDepolarization,radius);
		//basic.crossingTestAndVehicle(useDepolarization, radius);
		//basic.crossingTestMulti(false,subSteps,useDepolarization,radius);

//Curvature tests
		//CurvatureTests ct(sceneManager, radius);
		//ct.cylinderTest();
		//ct.symmetricDivergenceTest();
		
//Tunnel tests
		//CubeTunnelTests ctt(sceneManager,radius,useDepolarization);
		//ctt.cubeTunnelRDNIsotropic(100);
		//ctt.cubeTunnelRDN(100);
		//ctt.cubeTunnel(0);//No random
		//ctt.cubeTunnelSingleReceiver(0);//No random
		//ctt.cubeTunnelWithCurvedSimulation(0,1);

//RDN tests	
		
		//RDNTests rdn(sceneManager,radius,useDepolarization);
		//rdn.runDidascalouDielectricRDN(1);
		//rdn.runDidascalouConductorRDN(1);
		//rdn.runDidascalouDielectric(0);
		//rdn.runDidascalouDielectricMultipleReceivers(0);
		//rdn.runTest(test);
		//rdn.rdnPlaneTestMultiReceiver();
		//rdn.rdnPlaneTest();
		//rdn.freeSpace();
//Dudley tests
		//DudleyTests dud(sceneManager,radius,useDepolarization);
		//dud.runTest(test);		

//Anvers test
		//AnversTests anvers(sceneManager,radius,useDepolarization);
		//anvers.runTestsJuly21(test);
		//anvers.runTestsDepolarizationLeanAntenna(test, false); //true with antenna gain
		//anvers.runTestsBetweenLeanAntenna(test, true ); //true with antenna gain
		//anvers.runTests6GHz(test, true  ); //true with antenna gain
		//anvers.runTests(test);
		//anvers.runRDNSingleRx(true);
		//anvers.runRDNStraightTunnel(true);
//Anvers Mimo test
		//AnversMimo anvers(sceneManager,radius,useDepolarization);
		//anvers.runTests(test, outputPath);
	//	anvers.runLille(test, outputPath);
		
//AnversDoa test
		//AnversDoaTests anvers(sceneManager,radius,useDepolarization);
		//anvers.runTests(test);
		//anvers.runTests6GHz(test, false  ); //true with antenna gain


//RouxAP tests
		//RouxAPTests roux(sceneManager,radius,useDepolarization);
                //roux.runRDNCube();
		////roux.runADSingleRay();
		//roux.runTests(test);

//Diffraction tests
		//DiffractionTests dif(sceneManager,radius, useDepolarization);
		//dif.semiplaneDiffraction();
		//dif.semiplaneTotal();
		//dif.runSWCorner(false);
		//dif.runSWCorner(true); //With multitransmitter
		//dif.runSWAcuteCorner();
		//dif.runCrossing();
		//dif.addCompoundDynamicMeshes();
//AntennaGain test
		AntennaGainTests ag(sceneManager, radius);
		//ag.freeSpace(useDepolarization);
		//ag.freeSpaceRDN();
		ag.testAntennaOrientation();

//Cartagena GSM tests
		//CartagenaGSM ct(sceneManager,radius,useDepolarization);
		//ct.runTests(test);

//Parking tests
		//Parking p(sceneManager,radius, useDepolarization);
		//p.runTests(test, outputPath);
		//p.test();
//LoRa test
		//Lora lora(sceneManager,radius, useDepolarization);
		//lora.runTests(test);
		//lora.runAll(true,test, outputPath);
		//lora.runAllMoving(true,test, outputPath);
//Lille test
		//Lille lille(sceneManager,radius, useDepolarization);
		//lille.runBasestations(test,outputPath);
		//lille.testAntennaOrientation();
		//lille.runWaz(test,outputPath);
		//lille.runTests(test, outputPath );
//Exposure tests
//		Exposure exposure(sceneManager,radius, useDepolarization);
	//	exposure.run5GBeamSelection(test,outputPath);
//		exposure.run5GOneOperator(test,outputPath);
		delete sceneManager;

		return 0;
	}
	catch (optix::Exception& e) {
		std::cout << "main error occurred with error code "
			<< e.getErrorCode() << " and message "
			<< e.getErrorString() << std::endl;

		return 1;
	}
	catch (opal::Exception& e) {
		std::cout << "main error occurred with  message "
			<< e.getErrorString()
			<< std::endl;

		return 2;
	}

}



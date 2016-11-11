/*!
 * file
 * \brief a track reconstruction algorithm tester
 * \author Wolfgang Waltenberger
 * \date Sun 23 Sep 2001 23:32:56 CEST
 *
 * */

// $Id: traxx.cpp,v 1.3 2001/11/21 08:53:24 wwalten Exp $

#include <stdio.h>
#include <signal.h>
#include <unistd.h>
#include <sys/timeb.h>
#include <stdlib.h>
#include <vector.h>
#include <iostream>
#include <unistd.h>
#include <time.h>
#include <sys/times.h>

#include <iostream.h>
#include <fstream>
#include <CLHEP/config/CLHEP.h>
#include "CLHEP/Random/Randomize.h"

#include "../interface/DetElement.h"
#include "../interface/Detector.h"
#include "../interface/Hit.h"
#include "../interface/Plane.h"
#include "../interface/Track.h"
#include "../interface/defines.h"

#ifdef GL
#include <Inventor/Xt/SoXt.h>
#include <Inventor/nodes/SoMaterial.h>
#include <Inventor/Xt/SoXtRenderArea.h>
#include <Inventor/Xt/viewers/SoXtExaminerViewer.h>
#include <Inventor/nodes/SoDirectionalLight.h>
#include <Inventor/nodes/SoPerspectiveCamera.h>
#include <Inventor/nodes/SoSphere.h>
#include <Inventor/nodes/SoTransform.h>
#include <Inventor/Xt/SoXtGLWidget.h>
#include <Inventor/nodes/SoSelection.h>
#include <vector.h>
#endif

#include <string>
#include "getopt.h" // solaris doesnt have getopt.

#define DEFAULT_NTRACKS 1

#ifndef VERBOSE
#define VERBOSE 5
#endif

using namespace std; // solaris needs that.

string _main_argv;
extern int errno;

clock_t Timing[5];
clock_t userTiming[5];
clock_t sysTiming[5];
#ifdef DEBUG
clock_t lastTiming[5];
#endif
struct tms mytms;

char verbosity=VERBOSE;

#ifdef DEBUG
char watch_names[5][20]={"Simulation","FullReco","Full Overhead",
  "Covariance Matrix","Reference Track"};
#endif

inline clock_t wholeTiming( int i )
{
  return userTiming[i]+sysTiming[i];
};

void WatchStart( int i )
{
  Timing[i]-=times(&mytms);
  userTiming[i]-=mytms.tms_utime;
  sysTiming[i]-=mytms.tms_stime;
  #ifdef DEBUG
  lastTiming[i]=wholeTiming(i);
  MESSAGE(52,"[watch]    Started watch "+any2str(i) + " (" + watch_names[i] +
      ").");
  #endif
};

void WatchStop( int i )
{
  Timing[i]+=times(&mytms);
  userTiming[i]+=mytms.tms_utime;
  sysTiming[i]+=mytms.tms_stime;
  #ifdef DEBUG
  MESSAGE(53,"[watch]    Stopped watch " + any2str(i) + " (" + watch_names[i] +
      ").");
  MESSAGE(53,"[watch]       `---> " + any2str(wholeTiming(i) -lastTiming[i]));
  #endif
};

#ifdef GL
void invokeRestart(Widget, XtPointer ClientData, 
    XmPushButtonCallbackStruct *)
{
  string rstrt=_main_argv+" &";
  system (rstrt.c_str());
  exit (0);
};

void invokeRestart( int sig )
{
  string rstrt=_main_argv+" &";
  system (rstrt.c_str());
  exit(0);
};
#endif



void usage (char *name) {
  cout << "Usage: " << name << " [options]";
  cout << endl 
     << endl << "Operation: " << endl
     << "    -R, --riemann <levels>:  Riemann fit." << endl
     << "    -G, --global <levels>:   Global fit. " << endl
     << "    -K, --kalman <levels>:   Kalman filter." << endl
     << "  * <level> is the MS level: " << endl
     << "          0: no covariance matrix" << endl
     << "          1: approximate covariance matrix." << endl
     << "          2: exact covariance matrix" << endl
     << "          3: full covariance matrix" << endl
     << "  * levels can be concatenated. example:" << endl
     << "          --riemann 123" << endl
     << "  * the number of iterations can be given for Riemann fits" << endl
     << "    in forward detectors ( see --iterations )" << endl
     << "  * default mode: --riemann 012 --global 0123 --kal 01" << endl
     << endl << "Options: " << endl
     << "    -n, --ntracks <num>:  Number of tracks to simulate "
     << "(default: " << DEFAULT_NTRACKS << ")" << endl
     << "    -f, --detfile <file>: read detector configuration from"
     << " <file>." << endl
     << "                          This also changes the init file to be used"
     << endl
     << "                          default is "ETCDIR"/detector.conf" << endl
     << "    -i, --initfile <file>: read initial values from"
     << " <file>." << endl
     << "                          default is "ETCDIR"/detector.init.conf" << endl
     << "    -s, --seed <num>:     set seed to num." << endl
     << "    -a, --scatter <0/1/2>: scatter? 0: no 1: randomly 2: "
        " deterministically." << endl
     << "                          defaults to 1 (randomly)" << endl
     << "    -m, --smear <0/1/2>:  smear? 0: no 1: randomly 2: "
     << "deterministically" << endl
     << "                          defaults to 1 (randomly)." << endl
     << "    -I, --iterations <n>: Number of iterations for the Riemann Fit"
     << endl
     << "                          in a forward type detector; defaults to 1."
     << endl
     << endl
     << "    -p, --persistent:     persistent storage of tracks? Needs a lot"
     << " of memory" << endl 
     << "                          per track, but allows for precise timing" 
     << endl
     << "                          (default: false)" 
     << endl 
     << endl 
     << endl << "Output: " << endl
     << "    -c, --cms[=<pre>]:    write statistics (CMS coords) into files"
     << endl
     << "                          optionally prefixed with <pre>" << endl
     << "                          prefix also changes default timefile." <<
     endl
     << "    -d, --delphi[=<pre>]: write statistics in Delphi into files" 
     << endl
     << "                          optionally prefixed with <pre>" << endl
     << "                          prefix also changes default timefile."
     << endl
     << "    -T, --timing:         show the elapsed CPU time for the "
     << "                          for the different code regions into file." 
     << endl
     << "    -t, --timing[<file>]: write the elapsed CPU time for the "
     << endl
     << "                          for the different code regions into file." 
     << endl
     << "    -O, --rawout:               Write raw data to stdout. " << endl
     << "    -o <file>, --rawout=<file>: Write raw data to <file>."
            << " (Overrides -R)" << endl
     << "    -w, --write:                Write data human readably to stdout"
     << endl
     #ifdef GL
     << "    -g, --gui:            OpenGL graphics as output." << endl
     #endif
     << endl << "Informative output:" << endl
     << "    -h, --help:     This help screen, then exit." << endl
     << "    -V, --version:  Prints version number, then exit." << endl
     << "    -v, --verbose=<n>:  Change verbosity level (default is " <<
             VERBOSE << ")." << endl
     << "                        This changes the behavior in debug mode only."
     << endl << endl << "Bugs, wishes:" << endl 
     << "     Bug reports, feature requests, patches, etc. mailto:" << endl
     << "     Wolfgang Waltenberger <wwalten@hephy.oeaw.ac.at>" << endl;
  exit (1);
}

// these are global variables I needed for threading.
// now they may still be useful as globals.
namespace nmain {
  string initfile=ETCDIR"/detector.init.conf";
  int smear=1;
  int scatter=1;
  int iterations=1;
  char statistics=0;
  Detector *barrel;
  bool riemann[4]={0,0,0,0}, kalman[4]={0,0,0,0}, global[4]={0,0,0,0};
  int n_fits=0; // number of fits
  bool writelong=0,writeraw=0;
  string fname=""; // -r <file>
  FILE *file[12]; // faster than *file[n_fits]?
  #ifdef GL
  bool glout=0;
  SoSelection *root;
  #endif
};

using namespace nmain;

void evalTrack()
{
  // evaluateTrack
  // We simulate a track:
  Track track ( barrel );
  WatchStart ( 0 );
//  Timer.start();
  StateVector tv=track.Simulate( scatter, smear, initfile );
  FullVector fv[ 12 ];          // faster than fv [ n_fits ].
  WatchStop ( 0 );
  WatchStart ( 1 );
  int l=0;
  for (int j=0; j<4; j++) {
    if (riemann[j]) fv[l++]=track.RiemannFit ( j, iterations ) ;
    if (global[j])  fv[l++]=track.GlobalFit  ( j ) ;
    if (kalman[j])  fv[l++]=track.KalmanFilter (j);
  };
  WatchStop ( 1 );
  // Some output:
  if (writelong) { 
    track.write(); 
    for (int j=0; j<n_fits; j++) {
      cout << endl;
      fv[j].writeCMS();
      fv[j].writeDelphi();
    };
  };
  if (writeraw) {
    track.write_raw((char *) fname.c_str());
    for (int j=0; j<n_fits; j++) {
      fv[j].write_raw( (char *) fname.c_str() );
    };
  };
  if (statistics) {
    string line;
    for (int j=0;j<n_fits;j++) {
      if (statistics==2) {
        line=fv[j].DelphiValues ( tv );
      } else {
        line=fv[j].CMSValues ( tv );
      };
      /*
      char line[80];
      snprintf(line,80,"%.7e %.7e %.7e %.7e %.7e\n",a0()-tv.a0(),
        ,z()-tv.z(), theta()-tv.z(), psi(), kappa_cms());
      fputs ( myline, file[j] );
      */
      fputs(line.c_str(),file[j]);
      /*
      char blabbel[80];
      snprintf(blabbel,80,"%f %f\n",track.hits[11].Phi(),
          track.hits[11].R());
      fputs(blabbel, file[j]);
      */
    };
  };
  #ifdef GL
  if (glout) {
    // we draw the hits and the reconstructed track.
    track.GLdraw ( );
    fv[0].GLinitMaterial ( root ); // switch to blue.
    for (int j=0; j<n_fits; j++) {
      fv[j].GLdraw ( root , barrel->arclength() );
    };
  };
  #endif
};

void *evalNTracks_transient ( int n )
{
  MESSAGE(8,"[main]     evaluation will be transient.");
  if ( n > 99 ) {
    /*
    cerr << "8[                                                  ]" << endl
       << C_UP << C_RIGHT << "7";
       */
    // progress bar
    int step=n/50;
    for (int j=0; j<50; j++) {
      for ( int i=0; i<step; i++) {
        evalTrack();
      };
      cerr << C_INVERTED;
    };
    cerr << endl;
  } else {
    for ( int i=0; i<n; i++) {
      evalTrack();
    };
  };
  return NULL;
};

void *evalNTracks_persistent ( int n )
{
  vector < Track * > track;
  track.reserve(n);
  vector < StateVector > tv;
  tv.reserve(n);
  // Simulate all at once. better for timing.
  int step=n/50;
  MESSAGE(8,"[main]     evaluation will be persistent.");
  WatchStart ( 0 );
  for ( int i=0; i<n; i++ ) {
    track[i]=new Track ( barrel );
    tv.push_back( track[i]->Simulate( scatter, smear, initfile ) );
    if (!fmod (i,step)) {
      if (i==0) {
        cerr << "Simulate:    x";
      } else {
        cerr << "x";
      };
    };
  };
  WatchStop ( 0 );
  cout << endl << "Reconstruct: ";
  FullVector fv[12];          // faster than fv [ n_fits ].
  WatchStart ( 1 );
  for ( int i=0; i<n; i++ ) {
    if (!fmod (i,step)) cerr << "x";
    int l=0;
    for (int j=0; j<4; j++) {
      if (riemann[j]) fv[l++]=track[i]->RiemannFit ( j, iterations ) ;
      if (global[j])  fv[l++]=track[i]->GlobalFit  ( j ) ;
      if (kalman[j])  fv[l++]=track[i]->KalmanFilter (j);
    };
    if (writelong) { 
      track[i]->write(); 
      for (int j=0; j<n_fits; j++) {
        cout << endl;
        fv[j].writeCMS();
        fv[j].writeDelphi();
      };
    };
    if (writeraw) {
      track[i]->write_raw((char *) fname.c_str());
      for (int j=0; j<n_fits; j++) {
        fv[j].write_raw( (char *) fname.c_str() );
      };
    };
    if (statistics) {
      string line;
      for (int j=0;j<n_fits;j++) {
        if (statistics==2) {
          line=fv[j].DelphiValues ( tv[i] );
        } else {
          line=fv[j].CMSValues ( tv[i] );
        };
        fputs(line.c_str(),file[j]);
      };
    };
    #ifdef GL
    if (glout) {
      // we draw the hits and the reconstructed track.
      track[i]->GLdraw ( );
      fv[0].GLinitMaterial ( root ); // switch to blue.
      for (int j=0; j<n_fits; j++) {
        fv[j].GLdraw ( root , barrel->arclength() );
      };
    };
    #endif
  };
  cout << endl;
  WatchStop ( 1 );
  return NULL;
};

int main (int argc, char **argv)
{
  bool persistent=false;
  bool timing=false;
  bool timingfile=false;
  string timefile="timing";
  Timing[0]=0;
  Timing[1]=0;
  Timing[2]=0;
  Timing[3]=0;
  Timing[4]=0;
  signed long int longseed=0;
  bool seedgiven=false;
  
  int optret=0, ntracks=DEFAULT_NTRACKS;
  for (int i=0; i<argc;i++) {
    _main_argv+=argv[i];
    _main_argv+="  ";
  };
  
  string conffile=ETCDIR"/detector.conf"; // default config file
  string statfile="lastrun"; // prefix for statistics output
  
  
  bool mode_given=false;
  // long_options
  int options_index = 0;
  static struct option long_options[]=
  {
    {"help",      0, 0, 'h'},
    {"version",   0, 0, 'V'},
    {"verbose",   0, 0, 'v'},
    {"gui",       0, 0, 'g'},
    {"detfile",   1, 0, 'f'},
    {"initfile",  1, 0, 'i'},
    {"write",     0, 0, 'w'},
    {"rawout",    2, 0, 'o'},
    {"ntracks",   1, 0, 'n'},
    {"riemann",   1, 0, 'R'},
    {"global",    1, 0, 'G'},
    {"kalman",    1, 0, 'K'},
    {"cms",       2, 0, 'c'},
    {"delphi",    2, 0, 'd'},
    {"timing",    2, 0, 't'},
    {"seed",      1, 0, 's'},
    {"scatter",   1, 0, 'a'},
    {"smear",     1, 0, 'm'},
    {"iterations",1, 0, 'I'},
    {"persistent",0, 0, 'p'},
    {0 ,          0, 0 , 0}
  };
  
  while ((optret=getopt_long(argc,argv,
          "wR:G:K:gv:VhdOo:r:f:i:n:tTcds:e:m:a:pI:",long_options, \
          &options_index ))!=-1) {
    switch (optret) {
      case 'V': // --version
        cout << argv[0] << 
           " version: " << VERSION << endl 
           << "Built: " << __DATE__ << " " 
           << __TIME__ << endl;
        exit (0);
      case 'h': // --help
        usage(argv[0]);
        exit (0);
      case 'v': // --verbose
        verbosity=(char) atoi(optarg);
        break;
      case 'f':  // --detfile
        conffile=optarg;
        initfile=optarg;
        initfile.replace(initfile.find("conf"),4,"init.conf");
        break;
      case 'i':  // --initfile
        initfile=optarg;
        break;
      case 'O': // --rawout
        writeraw=true;
        fname="";
        break;
      case 'p': // --persistent
        persistent=true;
        break;
      case 'g': // --gui
        #ifdef GL
        glout=true;
        #else
        WARNING("No GL support compiled in");
        #endif
        break;
      case 'o':
        writeraw=true;
        if (optarg) {
//          cerr << "optarg=" << optarg;
//          fname=(char *) malloc(sizeof(optarg)+1);
//          snprintf(fname, sizeof(fname),optarg);
          fname=optarg;
        } else {
          fname="";
        };
        break;
      case 't': // --timing
        timingfile=true;
        if (optarg && optarg != "" ) timefile=optarg;
      case 'T':
        timing=true;
        break;
      case 'c': // --cms
        if (optarg) {
          statfile=optarg;
          timefile=statfile;
        };
        statistics=1;
        break;
      case 'd': // --delphi
        if (optarg) {
          statfile=optarg;
          timefile=statfile;
        };
        statistics=2;
        break;
      case 'R': // --riemann
        for (int i=0;;i++)
        {
          if (optarg[i]=='\0') break;
          int oa=(int) optarg[i]-48;
          if (oa<4 && oa>=0) {
            riemann[oa]=true;
            mode_given=true;
          };
        };
        break;
      case 'G': // --global
        for (int i=0;;i++)
        {
          if (optarg[i]=='\0') break;
          int oa=(int) optarg[i]-48;
          if (oa<4 && oa>=0) {
            global[oa]=true;
            mode_given=true;
          };
        };
        break;
      case 'K': // --kalman
        for (int i=0;;i++)
        {
          if (optarg[i]=='\0') break;
          int oa=(int) optarg[i]-48;
          if (oa<4 && oa>=0) {
            kalman[oa]=true;
            mode_given=true;
          };
        };
        break;
      case 'n': // --ntracks
        ntracks=atoi(optarg);
        if (ntracks <=0 ) {
          cerr << "Error: number of simulated tracks must be a "
                "positive nonzero integer " << endl;
          exit (-1);
        };
        break;
      case 's': // --seed
        longseed=atoi(optarg);
        seedgiven=true;
        break;
      case 'a': // --scatter
        scatter= optarg[0]-48;
        break;
      case 'I': // --iterations
        iterations = optarg[0]-48;
        break;
      case 'm': // --smear
        smear  = optarg[0]-48;
        break;
      case 'w': // --out
        writelong=true;
        break;
      default:
      cout << "Type '" << argv[0] << 
              " --help' to see the available options." << endl;
      exit (-1);
    };
  };
  
  if (!mode_given) { // these are the defaults.
    memset(riemann,true,3); // --riemann 012
    memset(global,true,4);  // --global 0123
    memset(kalman,true,2);  // --kalman 01
  };

  timeb seed;
  // initialize randnm number generator
  ftime(&seed);
  if (!seedgiven) {
    // abs value, cause negative long ints are somehow read in badly from
    // atoi
    longseed = (long int) (abs ( 1000 * seed.time + seed.millitm));
    srand48( longseed );
  };
  MESSAGE(8,"[system]   initializing random seed. Seed=" + any2str(longseed));
  HepRandom::setTheSeed ( longseed );
  //                  And now for the _real_ stuff
  //                 ------------------------------
  // We need a detector:
  barrel = new Detector ( conffile );
  
  if (writeraw) barrel->write_raw((char *) fname.c_str());
  
  #ifdef GL
  SoDB::init();
  // We draw the detector.
  root = new SoSelection;
  
  if (glout) {
    barrel->GLprologue ( root, "TinyCMS" );
    barrel->GLdraw ();
  };
  #endif
  
  for (int j=0;j<4;j++) n_fits+=riemann[j]+global[j]+kalman[j];
  // create the statistics files.
  if (statistics) {
    int j=0;
    char a[1];
    for (int i=0;i<4;i++) {
      sprintf(a,"%d",i);
      if (riemann[i]) {
        timefile+=".riemann"+any2str(i);
        file[j++]=fopen(
          (statfile+(string) ".riemann"+ (string) a).c_str(),"w");
      };
      if (global[i]) {
        timefile+=".global"+any2str(i);
        file[j++]=fopen(
          (statfile+(string) ".global"+ (string) a).c_str(),"w");
      };
      if (kalman[i]) {
        timefile+=".kalman"+any2str(i);
          file[j++]=fopen(
          (statfile+(string) ".kalman"+ (string) a).c_str(),"w");
      };
    };
  };
  
  // loop over ntracks.
  // two modi: transient, or persistent.
  // this simulates a track, reconstructs it, deletes the track,
  // simulates again ...
  if (persistent) {
    evalNTracks_persistent ( ntracks );
  } else {
    evalNTracks_transient ( ntracks );
  };
  
  #ifdef GL
  // we call GLepilogue to conclude the actual drawing.
  signal(SIGHUP,invokeRestart);
  
//  for (int j=0; j<n_fits;j++) { fclose(file[j]); };
  
  if (glout) barrel->GLepilogue ();
  #endif
  if (timing) {
    HepDouble cps = sysconf(_SC_CLK_TCK);
    if (timingfile) {
      timefile+=".t";
        string st= any2str( wholeTiming(0) ) + " " +
        any2str( wholeTiming(1) ) + " " +
        any2str( wholeTiming(2) ) + " " +
        any2str( wholeTiming(3) ) + " " +
        any2str( wholeTiming(4) );
      ofstream fout ( timefile.c_str() );
      MESSAGE(3,"[system]   wrote timing to file '"+timefile+"'.");
      fout << st;
    } else {
    cout << "Clocks per sec: " << cps << endl;
    cout << "Time elapsed in code region 0 ('simulation'):        "
         << ( wholeTiming(0) / cps) << " secs." << endl
         << "Time elapsed in code region 1 (''): "
         << ( wholeTiming(1) / cps) << " secs." << endl
         << "Time elapsed in code region 2 (''):           "
         << ( wholeTiming(2) / cps) << " secs." << endl
         << "Time elapsed in code region 3 (''):           "
         << ( wholeTiming(3) / cps) << " secs." << endl
         << "Time elapsed in code region 4 (''):           "
         << ( wholeTiming(4) / cps) << " secs." << endl
         << "----------------------------------------------------"
         << endl
         << "    => " << 1000 * wholeTiming(0) / cps / ntracks
         << " + " <<     1000 * wholeTiming(1) / cps / ntracks
         << " + " <<     1000 * wholeTiming(2) / cps / ntracks
         << " = " << 1000 * ( wholeTiming(2) + wholeTiming(0)) /
         cps / ntracks << " millisecs / track. "
           << endl << endl;
    cout << "  Real Time [Ticks]: " << Timing[0] << " " 
         << Timing[1] << " " << Timing[2] << endl;
    cout << "  User Time [Ticks]: " << userTiming[0] << " " 
         << userTiming[1] << " " << userTiming[2] << endl;
    cout << "System Time [Ticks]: " << sysTiming[0] << " " 
         << sysTiming[1] << " " << sysTiming[2] << endl;
    };
  };
  FINISH_COUTS();
};

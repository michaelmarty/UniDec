// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(Rcpp)]]
// [[Rcpp::plugins(unideclib)]]
#include <Rcpp.h>

#ifdef _WIN32
#include <windows.h>
#define DLL_HANDLE HMODULE
#define LOAD_DLL(name) LoadLibraryA(name)
#define GET_FN(handle, name) GetProcAddress(handle, name)
#else
#include <dlfcn.h>
#define DLL_HANDLE void*
#define LOAD_DLL(name) dlopen(name, RTLD_LAZY)
#define GET_FN(handle, name) dlsym(handle, name)
#endif

using namespace Rcpp;

#include "udstruct.h" // Assuming udstruct.h defines Config, Input, and function prototypes
#include "udinterface.hpp" // Exposes the InputVec and DeconVec classes in CPP

// [[Rcpp::export]]
SEXP create_config() {
    Rcpp::XPtr<Config> config(new Config(), true); // Declare local variable
    SetDefaultConfig(config); // Initialize with default values
    return config;
}

// [[Rcpp::export]]
SEXP create_inputvec() {
    return Rcpp::XPtr<InputVec>(new InputVec(), true);
}

// Function to set input data
// [[Rcpp::export]]
void set_input_data(Rcpp::XPtr<InputVec> input, Rcpp::XPtr<Config> config, NumericVector dataMZ, NumericVector dataInt) {
    input->dataMZ = std::vector<float>(dataMZ.begin(), dataMZ.end());
    input->dataInt = std::vector<float>(dataInt.begin(), dataInt.end());
    // Set lengthmz in config based on input data size
    config->lengthmz = input->dataMZ.size();
    // print lengthmz to verify
    printf("Input data set. lengthmz: %d\n", config->lengthmz);
}

// Function to set default config values
// [[Rcpp::export]]
void set_default_config(Rcpp::XPtr<Config> config) {
    SetDefaultConfig(config);
    printf("Default configuration set.\n");
}


// [[Rcpp::export]]
SEXP run_unidec(Rcpp::XPtr<Config> config, Rcpp::XPtr<InputVec> input) {
    DeconVec result = RunDeconvolution(config, input, 1);

    // Create S4 object of class "DeconResult"
    S4 s4res("DeconResult");
    s4res.slot("mass") = wrap(result.massaxis);
    s4res.slot("intensity") = wrap(result.massaxisval);
    // Set other slots as needed

    return s4res;
}

//None of this seemed to work, so going with XPtr approach
//// ----------- Expose Structs to R via Rcpp Module -----------
//
//RCPP_MODULE(unidec_mod) {
//    class_<Config>("Config")
//        .constructor()
//        // --- all scalar numeric fields, bind as .field() ---
//        .field("numit",              &Config::numit)
//        .field("numz",               &Config::numz)
//        .field("endz",               &Config::endz)
//        .field("startz",             &Config::startz)
//        .field("zsig",               &Config::zsig)
//        .field("psig",               &Config::psig)
//        .field("beta",               &Config::beta)
//        .field("mzsig",              &Config::mzsig)
//        .field("msig",               &Config::msig)
//        .field("molig",              &Config::molig)
//        .field("massub",             &Config::massub)
//        .field("masslb",             &Config::masslb)
//        .field("psfun",              &Config::psfun)
//        .field("zpsfun",             &Config::zpsfun)
//        .field("psmzthresh",         &Config::psmzthresh)
//        .field("mtabsig",            &Config::mtabsig)
//        .field("mflag",              &Config::mflag)
//        .field("massbins",           &Config::massbins)
//        .field("limitflag",          &Config::limitflag)
//        .field("psthresh",           &Config::psthresh)
//        .field("speedyflag",         &Config::speedyflag)
//        .field("linflag",            &Config::linflag)
//        .field("aggressiveflag",     &Config::aggressiveflag)
//        .field("adductmass",         &Config::adductmass)
//        .field("rawflag",            &Config::rawflag)
//        .field("nativezub",          &Config::nativezub)
//        .field("nativezlb",          &Config::nativezlb)
//        .field("poolflag",           &Config::poolflag)
//        .field("manualflag",         &Config::manualflag)
//        .field("intthresh",          &Config::intthresh)
//        .field("peakshapeinflate",   &Config::peakshapeinflate)
//        .field("fixedmassaxis",      &Config::fixedmassaxis)
//        .field("isotopemode",        &Config::isotopemode)
//        .field("filetype",           &Config::filetype)
//        .field("imflag",             &Config::imflag)
//
//        // IM Parameters
//        .field("dtsig",              &Config::dtsig)
//        .field("csig",               &Config::csig)
//        .field("ccsub",              &Config::ccsub)
//        .field("ccslb",              &Config::ccslb)
//        .field("ccsbins",            &Config::ccsbins)
//        .field("temp",               &Config::temp)
//        .field("press",              &Config::press)
//        .field("volt",               &Config::volt)
//        .field("tcal1",              &Config::tcal1)
//        .field("tcal2",              &Config::tcal2)
//        .field("tcal3",              &Config::tcal3)
//        .field("tcal4",              &Config::tcal4)
//        .field("twaveflag",          &Config::twaveflag)
//        .field("hmass",              &Config::hmass)
//        .field("to",                 &Config::to)
//        .field("len",                &Config::len)
//        .field("edc",                &Config::edc)
//        .field("nativeccsub",        &Config::nativeccsub)
//        .field("nativeccslb",        &Config::nativeccslb)
//        .field("baselineflag",       &Config::baselineflag)
//        .field("noiseflag",          &Config::noiseflag)
//        .field("zout",               &Config::zout)
//        .field("metamode",           &Config::metamode)
//        .field("minmz",              &Config::minmz)
//        .field("maxmz",              &Config::maxmz)
//        .field("mzres",              &Config::mzres)
//        .field("mzbins",             &Config::mzbins)
//        .field("bsub",               &Config::bsub)
//        .field("datareduction",      &Config::datareduction)
//        .field("peakwin",            &Config::peakwin)
//        .field("peakthresh",         &Config::peakthresh)
//        .field("exwindow",           &Config::exwindow)
//        .field("exchoice",           &Config::exchoice)
//        .field("exchoicez",          &Config::exchoicez)
//        .field("exthresh",           &Config::exthresh)
//        .field("exnorm",             &Config::exnorm)
//        .field("exnormz",            &Config::exnormz)
//        .field("peaknorm",           &Config::peaknorm)
//        .field("orbimode",           &Config::orbimode)
//        .field("datanorm",           &Config::datanorm)
//
//        // Experimental Parameters
//        .field("filterwidth",        &Config::filterwidth)
//        .field("zerolog",            &Config::zerolog)
//        .field("lengthmz",           &Config::lengthmz)
//        .field("mfilelen",           &Config::mfilelen)
//        .field("isolength",          &Config::isolength)
//
//        // DoubleDec Parameters
//        .field("doubledec",          &Config::doubledec)
//        .field("file_id",            &Config::file_id)
//        .field("silent",             &Config::silent)
//        .field("cdmsflag",           &Config::cdmsflag)
//
//        // Include set default parameters method from SetDefaultConfig(Config *config);
//        .method("SetDefaultConfig", SetDefaultConfig)
//    ;
//
//    // --- Input binding ---
//    class_<InputVec>("InputVec")
//        .constructor()
//        .field("dataMZ", &InputVec::dataMZ)
//        .field("dataInt", &InputVec::dataInt)
//        .field("testmasses", &InputVec::testmasses)
//        .field("nztab", &InputVec::nztab)
//        .field("mtab", &InputVec::mtab)
//        .field("barr", &InputVec::barr)
//        // Add custom methods here if needed
////        .method("setInput", &InputVec::setInput)
//        ;
//
//    // --- DeconVec binding ---
//    class_<DeconVec>("DeconVec")
//        .constructor()
//        .field("fitdat",     &DeconVec::fitdat)
//        .field("baseline",   &DeconVec::baseline)
//        .field("noise",      &DeconVec::noise)
//        .field("massgrid",   &DeconVec::massgrid)
//        .field("massaxis",   &DeconVec::massaxis)
//        .field("massaxisval",&DeconVec::massaxisval)
//        .field("blur",       &DeconVec::blur)
//        .field("newblur",    &DeconVec::newblur)
//        .field("peakx",      &DeconVec::peakx)
//        .field("peaky",      &DeconVec::peaky)
//        .field("dscores",    &DeconVec::dscores)
//        .field("starttab",   &DeconVec::starttab)
//        .field("endtab",     &DeconVec::endtab)
//        .field("mzdist",     &DeconVec::mzdist)
//        .field("rmzdist",    &DeconVec::rmzdist)
//        .field("error",      &DeconVec::error)
//        .field("rsquared",   &DeconVec::rsquared)
//        .field("iterations", &DeconVec::iterations)
//        .field("uniscore",   &DeconVec::uniscore)
//        .field("conv",       &DeconVec::conv)
//        .field("threshold",  &DeconVec::threshold)
//        .field("mlen",       &DeconVec::mlen)
//        .field("plen",       &DeconVec::plen)
//        .field("scanindex",  &DeconVec::scanindex)
//        ;
//
////    function("run_unidec", static_cast<int(*)(Config*, InputVec*)>(&run_unidec));
//    //function("run_unidec", &run_unidec);
//}








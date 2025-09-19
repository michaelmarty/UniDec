// [[Rcpp::plugins(cpp11)]]
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

// ----------- Expose Structs to R via Rcpp Module -----------

RCPP_MODULE(unidec_mod) {
    class_<Config>("Config")
        .constructor()
        // --- all scalar numeric fields, bind as .field() ---
        .field("numit",              &Config::numit)
        .field("numz",               &Config::numz)
        .field("endz",               &Config::endz)
        .field("startz",             &Config::startz)
        .field("zsig",               &Config::zsig)
        .field("psig",               &Config::psig)
        .field("beta",               &Config::beta)
        .field("mzsig",              &Config::mzsig)
        .field("msig",               &Config::msig)
        .field("molig",              &Config::molig)
        .field("massub",             &Config::massub)
        .field("masslb",             &Config::masslb)
        .field("psfun",              &Config::psfun)
        .field("zpsfun",             &Config::zpsfun)
        .field("psmzthresh",         &Config::psmzthresh)
        .field("mtabsig",            &Config::mtabsig)
        .field("mflag",              &Config::mflag)
        .field("massbins",           &Config::massbins)
        .field("limitflag",          &Config::limitflag)
        .field("psthresh",           &Config::psthresh)
        .field("speedyflag",         &Config::speedyflag)
        .field("linflag",            &Config::linflag)
        .field("aggressiveflag",     &Config::aggressiveflag)
        .field("adductmass",         &Config::adductmass)
        .field("rawflag",            &Config::rawflag)
        .field("nativezub",          &Config::nativezub)
        .field("nativezlb",          &Config::nativezlb)
        .field("poolflag",           &Config::poolflag)
        .field("manualflag",         &Config::manualflag)
        .field("intthresh",          &Config::intthresh)
        .field("peakshapeinflate",   &Config::peakshapeinflate)
        .field("fixedmassaxis",      &Config::fixedmassaxis)
        .field("isotopemode",        &Config::isotopemode)
        .field("filetype",           &Config::filetype)
        .field("imflag",             &Config::imflag)

        // IM Parameters
        .field("dtsig",              &Config::dtsig)
        .field("csig",               &Config::csig)
        .field("ccsub",              &Config::ccsub)
        .field("ccslb",              &Config::ccslb)
        .field("ccsbins",            &Config::ccsbins)
        .field("temp",               &Config::temp)
        .field("press",              &Config::press)
        .field("volt",               &Config::volt)
        .field("tcal1",              &Config::tcal1)
        .field("tcal2",              &Config::tcal2)
        .field("tcal3",              &Config::tcal3)
        .field("tcal4",              &Config::tcal4)
        .field("twaveflag",          &Config::twaveflag)
        .field("hmass",              &Config::hmass)
        .field("to",                 &Config::to)
        .field("len",                &Config::len)
        .field("edc",                &Config::edc)
        .field("nativeccsub",        &Config::nativeccsub)
        .field("nativeccslb",        &Config::nativeccslb)
        .field("baselineflag",       &Config::baselineflag)
        .field("noiseflag",          &Config::noiseflag)
        .field("zout",               &Config::zout)
        .field("metamode",           &Config::metamode)
        .field("minmz",              &Config::minmz)
        .field("maxmz",              &Config::maxmz)
        .field("mzres",              &Config::mzres)
        .field("mzbins",             &Config::mzbins)
        .field("bsub",               &Config::bsub)
        .field("datareduction",      &Config::datareduction)
        .field("peakwin",            &Config::peakwin)
        .field("peakthresh",         &Config::peakthresh)
        .field("exwindow",           &Config::exwindow)
        .field("exchoice",           &Config::exchoice)
        .field("exchoicez",          &Config::exchoicez)
        .field("exthresh",           &Config::exthresh)
        .field("exnorm",             &Config::exnorm)
        .field("exnormz",            &Config::exnormz)
        .field("peaknorm",           &Config::peaknorm)
        .field("orbimode",           &Config::orbimode)
        .field("datanorm",           &Config::datanorm)

        // Experimental Parameters
        .field("filterwidth",        &Config::filterwidth)
        .field("zerolog",            &Config::zerolog)
        .field("lengthmz",           &Config::lengthmz)
        .field("mfilelen",           &Config::mfilelen)
        .field("isolength",          &Config::isolength)

        // DoubleDec Parameters
        .field("doubledec",          &Config::doubledec)
        .field("file_id",            &Config::file_id)
        .field("silent",             &Config::silent)
        .field("cdmsflag",           &Config::cdmsflag)
    ;

    // --- Input binding ---
    class_<Input>("Input")
        .constructor()
        .field("dataMZ", &Input::dataMZ)
        .field("dataInt", &Input::dataInt)
        .field("testmasses", &Input::testmasses)
        .field("nztab", &Input::nztab)
        .field("mtab", &Input::mtab)
        .field("barr", &Input::barr)
        // Add custom methods here if needed
        ;

    // --- Decon binding ---
    class_<Decon>("Decon")
        .constructor()
        .field("fitdat",     &Decon::fitdat)
        .field("baseline",   &Decon::baseline)
        .field("noise",      &Decon::noise)
        .field("massgrid",   &Decon::massgrid)
        .field("massaxis",   &Decon::massaxis)
        .field("massaxisval",&Decon::massaxisval)
        .field("blur",       &Decon::blur)
        .field("newblur",    &Decon::newblur)
        .field("peakx",      &Decon::peakx)
        .field("peaky",      &Decon::peaky)
        .field("dscores",    &Decon::dscores)
        .field("starttab",   &Decon::starttab)
        .field("endtab",     &Decon::endtab)
        .field("mzdist",     &Decon::mzdist)
        .field("rmzdist",    &Decon::rmzdist)
        .field("error",      &Decon::error)
        .field("rsquared",   &Decon::rsquared)
        .field("iterations", &Decon::iterations)
        .field("uniscore",   &Decon::uniscore)
        .field("conv",       &Decon::conv)
        .field("threshold",  &Decon::threshold)
        .field("mlen",       &Decon::mlen)
        .field("plen",       &Decon::plen)
        .field("scanindex",  &Decon::scanindex)
        ;
}


// ----------- DLL Function Typedef for Example -----------

typedef int (*run_unidec_core_t)(Config*, Input*, int, int);

// ----------- Example Function to Call DLL Function -----------

// [[Rcpp::export]]
int call_run_unidec_core(std::string dll_path, SEXP configR, SEXP inputR) {
    // Retrieve pointers to existing Config and Input objects
    Rcpp::XPtr<Config> config(configR);
    Rcpp::XPtr<Input> input(inputR);

    // Load DLL
    DLL_HANDLE dh = LOAD_DLL(dll_path.c_str());
    if (!dh)
        Rcpp::stop("Failed to load DLL: %s", dll_path);

    // Get function pointer
    run_unidec_core_t fn = (run_unidec_core_t)GET_FN(dh, "run_unidec_core");
    if (!fn)
        Rcpp::stop("Failed to locate run_unidec_core in DLL");

    // Call function
    int result = fn(config.get(), input.get(), 0, 0);

    // Unload DLL
#ifdef _WIN32
    FreeLibrary(dh);
#else
    dlclose(dh);
#endif

    return result;
}

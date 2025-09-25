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

//void S4_to_Config(SEXP configS4, Config &config) {
//    config.numit              = Rcpp::as<int>(configS4.slot("numit"));
//    config.numz               = Rcpp::as<int>(configS4.slot("numz"));
//    config.endz               = Rcpp::as<int>(configS4.slot("endz"));
//    config.startz             = Rcpp::as<int>(configS4.slot("startz"));
//    config.zsig               = Rcpp::as<float>(configS4.slot("zsig"));
//    config.psig               = Rcpp::as<float>(configS4.slot("psig"));
//    config.beta               = Rcpp::as<float>(configS4.slot("beta"));
//    config.mzsig              = Rcpp::as<float>(configS4.slot("mzsig"));
//    config.msig               = Rcpp::as<float>(configS4.slot("msig"));
//    config.molig              = Rcpp::as<float>(configS4.slot("molig"));
//    config.massub             = Rcpp::as<float>(configS4.slot("massub"));
//    config.masslb             = Rcpp::as<float>(configS4.slot("masslb"));
//    config.psfun              = Rcpp::as<int>(configS4.slot("psfun"));
//    config.zpsfun             = Rcpp::as<int>(configS4.slot("zpsfun"));
//    config.psmzthresh         = Rcpp::as<float>(configS4.slot("psmzthresh"));
//    config.mtabsig            = Rcpp::as<float>(configS4.slot("mtabsig"));
//    config.mflag              = Rcpp::as<int>(configS4.slot("mflag"));
//    config.massbins           = Rcpp::as<float>(configS4.slot("massbins"));
//    config.limitflag          = Rcpp::as<int>(configS4.slot("limitflag"));
//    config.psthresh           = Rcpp::as<float>(configS4.slot("psthresh"));
//    config.speedyflag         = Rcpp::as<int>(configS4.slot("speedyflag"));
//    config.linflag            = Rcpp::as<int>(configS4.slot("linflag"));
//    config.aggressiveflag     = Rcpp::as<int>(configS4.slot("aggressiveflag"));
//    config.adductmass         = Rcpp::as<float>(configS4.slot("adductmass"));
//    config.rawflag            = Rcpp::as<int>(configS4.slot("rawflag"));
//    config.nativezub          = Rcpp::as<int>(configS4.slot("nativezub"));
//    config.nativezlb          = Rcpp::as<int>(configS4.slot("nativezlb"));
//    config.poolflag           = Rcpp::as<int>(configS4.slot("poolflag"));
//    config.manualflag         = Rcpp::as<int>(configS4.slot("manualflag"));
//    config.intthresh          = Rcpp::as<float>(configS4.slot("intthresh"));
//    config.peakshapeinflate   = Rcpp::as<float>(configS4.slot("peakshapeinflate"));
//    config.fixedmassaxis      = Rcpp::as<int>(configS4.slot("fixedmassaxis"));
//    config.isotopemode        = Rcpp::as<int>(configS4.slot("isotopemode"));
//    config.filetype           = Rcpp::as<int>(configS4.slot("filetype"));
//    config.imflag             = Rcpp::as<int>(configS4.slot("imflag"));
//    // IM Parameters
//    config.dtsig              = Rcpp::as<float>(configS4.slot("dtsig"));
//    config.csig               = Rcpp::as<float>(configS4.slot("csig"));
//    config.ccsub              = Rcpp::as<float>(configS4.slot("ccsub"));
//    config.ccslb              = Rcpp::as<float>(configS4.slot("ccslb"));
//    config.ccsbins            = Rcpp::as<float>(configS4.slot("ccsbins"));
//    config.temp               = Rcpp::as<float>(configS4.slot("temp"));
//    config.press              = Rcpp::as<float>(configS4.slot("press"));
//    config.volt               = Rcpp::as<float>(configS4.slot("volt"));
//    config.tcal1              = Rcpp::as<float>(configS4.slot("tcal1"));
//    config.tcal2              = Rcpp::as<float>(configS4.slot("tcal2"));
//    config.tcal3              = Rcpp::as<float>(configS4.slot("tcal3"));
//    config.tcal4              = Rcpp::as<float>(configS4.slot("tcal4"));
//    config.twaveflag          = Rcpp::as<int>(configS4.slot("twaveflag"));
//    config.hmass              = Rcpp::as<float>(configS4.slot("hmass"));
//    config.to                 = Rcpp::as<float>(configS4.slot("to"));
//    config.len                = Rcpp::as<float>(configS4.slot("len"));
//    config.edc                = Rcpp::as<float>(configS4.slot("edc"));
//    config.nativeccsub        = Rcpp::as<int>(configS4.slot("nativeccsub"));
//    config.nativeccslb        = Rcpp::as<int>(configS4.slot("nativeccslb"));
//    config.baselineflag       = Rcpp::as<int>(configS4.slot("baselineflag"));
//    config.noiseflag          = Rcpp::as<int>(configS4.slot("noiseflag"));
//    config.zout               = Rcpp::as<int>(configS4.slot("zout"));
//    config.metamode           = Rcpp::as<int>(configS4.slot("metamode"));
//    config.minmz              = Rcpp::as<float>(configS4.slot("minmz"));
//    config.maxmz              = Rcpp::as<float>(configS4.slot("maxmz"));
//    config.mzres              = Rcpp::as<float>(configS4.slot("mzres"));
//    config.mzbins             = Rcpp::as<float>(configS4.slot("mzbins"));
//    config.bsub               = Rcpp::as<float>(configS4.slot("bsub"));
//    config.datareduction      = Rcpp::as<int>(configS4.slot("datareduction"));
//    config.peakwin            = Rcpp::as<int>(configS4.slot("peakwin"));
//    config.peakthresh         = Rcpp::as<float>(configS4.slot("peakthresh"));
//    config.exwindow           = Rcpp::as<int>(configS4.slot("exwindow"));
//    config.exchoice           = Rcpp::as<int>(configS4.slot("exchoice"));
//    config.exchoicez          = Rcpp::as<int>(configS4.slot("exchoicez"));
//    config.exthresh           = Rcpp::as<float>(configS4.slot("exthresh"));
//    config.exnorm             = Rcpp::as<int>(configS4.slot("exnorm"));
//    config.exnormz            = Rcpp::as<int>(configS4.slot("exnormz"));
//    config.peaknorm           = Rcpp::as<int>(configS4.slot("peaknorm"));
//    config.orbimode           = Rcpp::as<int>(configS4.slot("orbimode"));
//    config.datanorm           = Rcpp::as<int>(configS4.slot("datanorm"));
//    // Experimental Parameters
//    config.filterwidth        = Rcpp::as<int>(configS4.slot("filterwidth"));
//    config.zerolog            = Rcpp::as<int>(configS4.slot("zerolog"));
//    config.lengthmz           = Rcpp::as<int>(configS4.slot("lengthmz"));
//    config.mfilelen           = Rcpp::as<int>(configS4.slot("mfilelen"));
//    config.isolength          = Rcpp::as<int>(configS4.slot("isolength"));
//    // DoubleDec Parameters
//    config.doubledec          = Rcpp::as<int>(configS4.slot("doubledec"));
//    config.file_id            = Rcpp::as<int>(configS4.slot("file_id"));
//    config.silent             = Rcpp::as<int>(configS4.slot("silent"));
//    config.cdmsflag           = Rcpp::as<int>(configS4.slot("cdmsflag"));
//}

#include <Rcpp.h>
using namespace Rcpp;

#include <Rcpp.h>
using namespace Rcpp;

void SEXP_to_Config(SEXP configSEXP, Config &config) {
    Rcpp::List configS4(configSEXP); // treat SEXP as a named list

    config.numit              = Rcpp::as<int>(configS4["numit"]);
    config.numz               = Rcpp::as<int>(configS4["numz"]);
    config.endz               = Rcpp::as<int>(configS4["endz"]);
    config.startz             = Rcpp::as<int>(configS4["startz"]);
    config.zsig               = Rcpp::as<float>(configS4["zsig"]);
    config.psig               = Rcpp::as<float>(configS4["psig"]);
    config.beta               = Rcpp::as<float>(configS4["beta"]);
    config.mzsig              = Rcpp::as<float>(configS4["mzsig"]);
    config.msig               = Rcpp::as<float>(configS4["msig"]);
    config.molig              = Rcpp::as<float>(configS4["molig"]);
    config.massub             = Rcpp::as<float>(configS4["massub"]);
    config.masslb             = Rcpp::as<float>(configS4["masslb"]);
    config.psfun              = Rcpp::as<int>(configS4["psfun"]);
    config.zpsfun             = Rcpp::as<int>(configS4["zpsfun"]);
    config.psmzthresh         = Rcpp::as<float>(configS4["psmzthresh"]);
    config.mtabsig            = Rcpp::as<float>(configS4["mtabsig"]);
    config.mflag              = Rcpp::as<int>(configS4["mflag"]);
    config.massbins           = Rcpp::as<float>(configS4["massbins"]);
    config.limitflag          = Rcpp::as<int>(configS4["limitflag"]);
    config.psthresh           = Rcpp::as<float>(configS4["psthresh"]);
    config.speedyflag         = Rcpp::as<int>(configS4["speedyflag"]);
    config.linflag            = Rcpp::as<int>(configS4["linflag"]);
    config.aggressiveflag     = Rcpp::as<int>(configS4["aggressiveflag"]);
    config.adductmass         = Rcpp::as<float>(configS4["adductmass"]);
    config.rawflag            = Rcpp::as<int>(configS4["rawflag"]);
    config.nativezub          = Rcpp::as<int>(configS4["nativezub"]);
    config.nativezlb          = Rcpp::as<int>(configS4["nativezlb"]);
    config.poolflag           = Rcpp::as<int>(configS4["poolflag"]);
    config.manualflag         = Rcpp::as<int>(configS4["manualflag"]);
    config.intthresh          = Rcpp::as<float>(configS4["intthresh"]);
    config.peakshapeinflate   = Rcpp::as<float>(configS4["peakshapeinflate"]);
    config.fixedmassaxis      = Rcpp::as<int>(configS4["fixedmassaxis"]);
    config.isotopemode        = Rcpp::as<int>(configS4["isotopemode"]);
    config.filetype           = Rcpp::as<int>(configS4["filetype"]);
    config.imflag             = Rcpp::as<int>(configS4["imflag"]);

    // IM Parameters
    config.dtsig              = Rcpp::as<float>(configS4["dtsig"]);
    config.csig               = Rcpp::as<float>(configS4["csig"]);
    config.ccsub              = Rcpp::as<float>(configS4["ccsub"]);
    config.ccslb              = Rcpp::as<float>(configS4["ccslb"]);
    config.ccsbins            = Rcpp::as<float>(configS4["ccsbins"]);
    config.temp               = Rcpp::as<float>(configS4["temp"]);
    config.press              = Rcpp::as<float>(configS4["press"]);
    config.volt               = Rcpp::as<float>(configS4["volt"]);
    config.tcal1              = Rcpp::as<float>(configS4["tcal1"]);
    config.tcal2              = Rcpp::as<float>(configS4["tcal2"]);
    config.tcal3              = Rcpp::as<float>(configS4["tcal3"]);
    config.tcal4              = Rcpp::as<float>(configS4["tcal4"]);
    config.twaveflag          = Rcpp::as<int>(configS4["twaveflag"]);
    config.hmass              = Rcpp::as<float>(configS4["hmass"]);
    config.to                 = Rcpp::as<float>(configS4["to"]);
    config.len                = Rcpp::as<float>(configS4["len"]);
    config.edc                = Rcpp::as<float>(configS4["edc"]);
    config.nativeccsub        = Rcpp::as<int>(configS4["nativeccsub"]);
    config.nativeccslb        = Rcpp::as<int>(configS4["nativeccslb"]);
    config.baselineflag       = Rcpp::as<int>(configS4["baselineflag"]);
    config.noiseflag          = Rcpp::as<int>(configS4["noiseflag"]);
    config.zout               = Rcpp::as<int>(configS4["zout"]);
    config.metamode           = Rcpp::as<int>(configS4["metamode"]);
    config.minmz              = Rcpp::as<float>(configS4["minmz"]);
    config.maxmz              = Rcpp::as<float>(configS4["maxmz"]);
    config.mzres              = Rcpp::as<float>(configS4["mzres"]);
    config.mzbins             = Rcpp::as<float>(configS4["mzbins"]);
    config.bsub               = Rcpp::as<float>(configS4["bsub"]);
    config.datareduction      = Rcpp::as<int>(configS4["datareduction"]);
    config.peakwin            = Rcpp::as<int>(configS4["peakwin"]);
    config.peakthresh         = Rcpp::as<float>(configS4["peakthresh"]);
    config.exwindow           = Rcpp::as<int>(configS4["exwindow"]);
    config.exchoice           = Rcpp::as<int>(configS4["exchoice"]);
    config.exchoicez          = Rcpp::as<int>(configS4["exchoicez"]);
    config.exthresh           = Rcpp::as<float>(configS4["exthresh"]);
    config.exnorm             = Rcpp::as<int>(configS4["exnorm"]);
    config.exnormz            = Rcpp::as<int>(configS4["exnormz"]);
    config.peaknorm           = Rcpp::as<int>(configS4["peaknorm"]);
    config.orbimode           = Rcpp::as<int>(configS4["orbimode"]);
    config.datanorm           = Rcpp::as<int>(configS4["datanorm"]);

    // Experimental Parameters
    config.filterwidth        = Rcpp::as<int>(configS4["filterwidth"]);
    config.zerolog            = Rcpp::as<int>(configS4["zerolog"]);
    config.lengthmz           = Rcpp::as<int>(configS4["lengthmz"]);
    config.mfilelen           = Rcpp::as<int>(configS4["mfilelen"]);
    config.isolength          = Rcpp::as<int>(configS4["isolength"]);

    // DoubleDec Parameters
    config.doubledec          = Rcpp::as<int>(configS4["doubledec"]);
    config.file_id            = Rcpp::as<int>(configS4["file_id"]);
    config.silent             = Rcpp::as<int>(configS4["silent"]);
    config.cdmsflag           = Rcpp::as<int>(configS4["cdmsflag"]);
}


//// S4-based function, NOT exported with sourceCpp
Rcpp::S4 run_unidec(SEXP configS4, SEXP inputS4) {
    // Create a blank InputVec and Config object
    InputVec input = InputVec();
    Config config = Config();

    // Convert S4 to Config and InputVec
    SEXP_to_Config(configS4, config);

//    Rcpp::XPtr<Config> config(configS4);
//    Rcpp::XPtr<InputVec> input(inputS4);
    DeconVec result = RunDeconvolution(config, input, 0, 0);
    // Convrt DeconVec to S4
    Rcpp::S4 resultS4("DeconVec");
    resultS4.slot("fitdat")      = Rcpp::NumericVector(result.fitdat.begin(), result.fitdat.end());
    resultS4.slot("baseline")    = Rcpp::NumericVector(result.baseline.begin(), result.baseline.end());
    resultS4.slot("noise")      = Rcpp::NumericVector(result.noise.begin(), result.noise.end());
    resultS4.slot("massgrid")   = Rcpp::NumericVector(result.massgrid.begin(), result.massgrid.end());
    resultS4.slot("massaxis")   = Rcpp::NumericVector(result.massaxis.begin(), result.massaxis.end());
    resultS4.slot("massaxisval")= Rcpp::NumericVector(result.massaxisval.begin(), result.massaxisval.end());
    resultS4.slot("blur")       = Rcpp::NumericVector(result.blur.begin(), result.blur.end());
    resultS4.slot("newblur")    = Rcpp::NumericVector(result.newblur.begin(), result.newblur.end());
    resultS4.slot("peakx")      = Rcpp::NumericVector(result.peakx.begin(), result.peakx.end());
    resultS4.slot("peaky")      = Rcpp::NumericVector(result.peaky.begin(), result.peaky.end());
    resultS4.slot("dscores")    = Rcpp::NumericVector(result.dscores.begin(), result.dscores.end());
    resultS4.slot("starttab")   = Rcpp::IntegerVector(result.starttab.begin(), result.starttab.end());
    resultS4.slot("endtab")     = Rcpp::IntegerVector(result.endtab.begin(), result.endtab.end());
    resultS4.slot("mzdist")     = Rcpp::NumericVector(result.mzdist.begin(), result.mzdist.end());
    resultS4.slot("rmzdist")    = Rcpp::NumericVector(result.rmzdist.begin(), result.rmzdist.end());
    resultS4.slot("error")      = result.error;
    resultS4.slot("rsquared")   = result.rsquared;
    resultS4.slot("iterations") = result.iterations;
    resultS4.slot("uniscore")   = result.uniscore;
    resultS4.slot("conv")       = result.conv;
    resultS4.slot("threshold")  = result.threshold;
    resultS4.slot("mlen")       = result.mlen;
    resultS4.slot("plen")       = result.plen;
    resultS4.slot("scanindex")  = result.scanindex;

    return resultS4;
}
//
//// // [[Rcpp::export]]
//int run_unidec(Rcpp::SEXP configR, Rcpp::SEXP inputR) {
//    Rcpp::XPtr<Config> config(configR);
//    Rcpp::XPtr<InputVec> input(inputR);
//
//    // Pass references, not pointers
//    DeconVec result = RunDeconvolution(*config, *input, 0, 0);
//
////    return Rcpp::List::create(
////        Rcpp::Named("massaxis") = Rcpp::NumericVector(result.massaxis.begin(), result.massaxis.end()),
////        Rcpp::Named("massaxisval") = Rcpp::NumericVector(result.massaxisval.begin(), result.massaxisval.end())
////        // Add other fields as needed
////    );
//return 0;
//}

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

        // Include set default parameters method from SetDefaultConfig(Config *config);
        .method("SetDefaultConfig", SetDefaultConfig)
    ;

    // --- Input binding ---
    class_<InputVec>("InputVec")
        .constructor()
        .field("dataMZ", &InputVec::dataMZ)
        .field("dataInt", &InputVec::dataInt)
        .field("testmasses", &InputVec::testmasses)
        .field("nztab", &InputVec::nztab)
        .field("mtab", &InputVec::mtab)
        .field("barr", &InputVec::barr)
        // Add custom methods here if needed
        .method("setInput", &InputVec::setInput)
        ;

    // --- DeconVec binding ---
    class_<DeconVec>("DeconVec")
        .constructor()
        .field("fitdat",     &DeconVec::fitdat)
        .field("baseline",   &DeconVec::baseline)
        .field("noise",      &DeconVec::noise)
        .field("massgrid",   &DeconVec::massgrid)
        .field("massaxis",   &DeconVec::massaxis)
        .field("massaxisval",&DeconVec::massaxisval)
        .field("blur",       &DeconVec::blur)
        .field("newblur",    &DeconVec::newblur)
        .field("peakx",      &DeconVec::peakx)
        .field("peaky",      &DeconVec::peaky)
        .field("dscores",    &DeconVec::dscores)
        .field("starttab",   &DeconVec::starttab)
        .field("endtab",     &DeconVec::endtab)
        .field("mzdist",     &DeconVec::mzdist)
        .field("rmzdist",    &DeconVec::rmzdist)
        .field("error",      &DeconVec::error)
        .field("rsquared",   &DeconVec::rsquared)
        .field("iterations", &DeconVec::iterations)
        .field("uniscore",   &DeconVec::uniscore)
        .field("conv",       &DeconVec::conv)
        .field("threshold",  &DeconVec::threshold)
        .field("mlen",       &DeconVec::mlen)
        .field("plen",       &DeconVec::plen)
        .field("scanindex",  &DeconVec::scanindex)
        ;

    function("run_unidec", static_cast<Rcpp::S4(*)(SEXP, SEXP)>(&run_unidec));
}








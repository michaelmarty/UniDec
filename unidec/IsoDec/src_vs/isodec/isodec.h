#pragma once

#ifdef ISODEC_EXPORTS
#define ISODEC_API __declspec(dllexport)
#else
#define ISODEC_API __declspec(dllimport)
#endif

extern "C" ISODEC_API void run();
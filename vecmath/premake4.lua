project "vecmath"
  kind "StaticLib"
  language "C++"
  files{
    "**.cpp",
    "**.h"
  }
  includedirs{"include"}
  configuration {"Debug"}
    defines{"Debug"}
    flags{"Symbols"}
  
  configuration {"Release"}
    defines{"NDebug"}
    flags{"Optimize"}
  

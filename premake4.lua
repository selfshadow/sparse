solution "sparse"
   configurations { "Debug", "Release" }

   project "ksvd_image"
      kind "ConsoleApp"
      language "C++"
      files { "**.h", "**.cpp", "**.c" }

      configuration "Debug"
         defines { "DEBUG" }
         flags { "Symbols" }

      configuration "Release"
         defines { "NDEBUG" }
         flags { "Optimize" }

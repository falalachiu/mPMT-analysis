/*
Author:    Mohit Gola 5th February 2024
*/

#include "constructInSitumPMT.hh"
#include "construction.hh"
#include "constructionPMT.hh"
#include "inSituPMT.hh"
#include "G4MultiUnion.hh"
#include "G4SubtractionSolid.hh"
#include "WCSimTuningParameters.hh"
#include "G4LogicalSkinSurface.hh"
#include <algorithm>
#include <vector>
#include <tuple>
ConstructInSitumPMT::ConstructInSitumPMT() {}

ConstructInSitumPMT::~ConstructInSitumPMT() {}

G4VPhysicalVolume *ConstructInSitumPMT::Construct()
{

  ///////////////////construct material/////////////////////

  G4NistManager *nist = G4NistManager::Instance();

  const G4int NUMENTRIES_water = 60;

  G4double ENERGY_water[NUMENTRIES_water] =
      {1.56962 * eV, 1.58974 * eV, 1.61039 * eV, 1.63157 * eV,
       1.65333 * eV, 1.67567 * eV, 1.69863 * eV, 1.72222 * eV,
       1.74647 * eV, 1.77142 * eV, 1.7971 * eV, 1.82352 * eV,
       1.85074 * eV, 1.87878 * eV, 1.90769 * eV, 1.93749 * eV,
       1.96825 * eV, 1.99999 * eV, 2.03278 * eV, 2.06666 * eV,
       2.10169 * eV, 2.13793 * eV, 2.17543 * eV, 2.21428 * eV,
       2.25454 * eV, 2.29629 * eV, 2.33962 * eV, 2.38461 * eV,
       2.43137 * eV, 2.47999 * eV, 2.53061 * eV, 2.58333 * eV,
       2.63829 * eV, 2.69565 * eV, 2.75555 * eV, 2.81817 * eV,
       2.88371 * eV, 2.95237 * eV, 3.02438 * eV, 3.09999 * eV,
       3.17948 * eV, 3.26315 * eV, 3.35134 * eV, 3.44444 * eV,
       3.54285 * eV, 3.64705 * eV, 3.75757 * eV, 3.87499 * eV,
       3.99999 * eV, 4.13332 * eV, 4.27585 * eV, 4.42856 * eV,
       4.59258 * eV, 4.76922 * eV, 4.95999 * eV, 5.16665 * eV,
       5.39129 * eV, 5.63635 * eV, 5.90475 * eV, 6.19998 * eV};

  G4double RINDEX_air[NUMENTRIES_water] =
      {1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
       1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
       1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
       1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
       1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
       1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
       1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
       1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
       1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
       1.0, 1.0, 1.0, 1.0, 1.0, 1.0};

  G4double ENERGY_SilGel[87] = {4.678776674 * eV, 4.592027283 * eV, 4.508478709 * eV, 4.427958491 * eV, 4.35030587 * eV,
                                4.275369924 * eV, 4.202452758 * eV, 4.13256199 * eV, 4.064996939 * eV, 3.999137915 * eV,
                                3.935906953 * eV, 3.87420652 * eV, 3.814907862 * eV, 3.756987416 * eV, 3.700835849 * eV,
                                3.646372728 * eV, 3.593524966 * eV, 3.542221472 * eV, 3.492396167 * eV, 3.443986692 * eV,
                                3.396571472 * eV, 3.350829881 * eV, 3.305993162 * eV, 3.262706518 * eV, 3.220245191 * eV,
                                3.178906699 * eV, 3.138646881 * eV, 3.099424803 * eV, 3.061201121 * eV, 3.023938677 * eV,
                                2.987602577 * eV, 2.95188653 * eV, 2.917309072 * eV, 2.883302061 * eV, 2.850107745 * eV,
                                2.817696981 * eV, 2.786043724 * eV, 2.7551212 * eV, 2.72490498 * eV, 2.695146764 * eV,
                                2.666278206 * eV, 2.637833975 * eV, 2.610016557 * eV, 2.583011752 * eV, 2.556384369 * eV,
                                2.53013002 * eV, 2.504627471 * eV, 2.479660085 * eV, 2.455025471 * eV, 2.430901116 * eV,
                                2.407447973 * eV, 2.384294297 * eV, 2.361606176 * eV, 2.339203259 * eV, 2.317409616 * eV,
                                2.295881811 * eV, 2.274931174 * eV, 2.254228848 * eV, 2.233923745 * eV, 2.213856184 * eV,
                                2.194315186 * eV, 2.175138494 * eV, 2.156177117 * eV, 2.137566187 * eV, 2.119296781 * eV,
                                2.10135941 * eV, 2.083745787 * eV, 2.06631904 * eV, 2.049203275 * eV, 2.032514943 * eV,
                                2.01587474 * eV, 1.999646316 * eV, 1.983698791 * eV, 1.967909817 * eV, 1.952391965 * eV,
                                1.937138174 * eV, 1.92214235 * eV, 1.90739769 * eV, 1.892792856 * eV, 1.878535228 * eV,
                                1.864408987 * eV, 1.850413879 * eV, 1.836747218 * eV, 1.823203874 * eV, 1.809879596 * eV,
                                1.796768999 * eV, 1.783867346 * eV};

  G4double ABSORPTION_SilGel[87] = {0.038149438 * m, 0.042804385 * m, 0.052016858 * m, 0.065315823 * m, 0.082271384 * m,
                                    0.102492373 * m, 0.125624001 * m, 0.151345579 * m, 0.179368306 * m, 0.209433124 * m,
                                    0.241308636 * m, 0.274789092 * m, 0.309692442 * m, 0.345858458 * m, 0.383146918 * m,
                                    0.421435858 * m, 0.460619889 * m, 0.500608586 * m, 0.541324938 * m, 0.582703865 * m,
                                    0.624690804 * m, 0.667240357 * m, 0.710315014 * m, 0.753883932 * m, 0.797921786 * m,
                                    0.842407687 * m, 0.887324167 * m, 0.932656224 * m, 0.978390441 * m, 1.024514168 * m,
                                    1.071014774 * m, 1.117878956 * m, 1.165092126 * m, 1.212637859 * m, 1.260497406 * m,
                                    1.308649275 * m, 1.357068882 * m, 1.405728262 * m, 1.45459585 * m, 1.503636328 * m,
                                    1.552810543 * m, 1.602075479 * m, 1.65138431 * m, 1.700686509 * m, 1.74992803 * m,
                                    1.799051552 * m, 1.847996793 * m, 1.896700886 * m, 1.945098826 * m, 1.993123983 * m,
                                    2.040708675 * m, 2.08778482 * m, 2.134284638 * m, 2.180141438 * m, 2.225290453 * m,
                                    2.269669756 * m, 2.313221236 * m, 2.355891637 * m, 2.397633674 * m, 2.438407204 * m,
                                    2.478180472 * m, 2.516931417 * m, 2.55464905 * m, 2.591334893 * m, 2.627004493 * m,
                                    2.661688988 * m, 2.695436756 * m, 2.72831512 * m, 2.760412122 * m, 2.791838362 * m,
                                    2.822728908 * m, 2.853245269 * m, 2.883577433 * m, 2.913945973 * m, 2.944604226 * m,
                                    2.975840522 * m, 3.0079805 * m, 3.041389473 * m, 3.076474872 * m, 3.113688745 * m,
                                    3.153530334 * m, 3.196548711 * m, 3.243345483 * m, 3.29457756 * m, 3.350959998 * m,
                                    3.413268897 * m, 3.482344375 * m};

  G4double RINDEX_SilGel[87] = {1.608013009, 1.591127095, 1.575602811, 1.561337877, 1.548237552,
                                1.536214093, 1.52518625, 1.515078789, 1.505822046, 1.497351509,
                                1.489607425, 1.482534437, 1.476081242, 1.470200274, 1.464847411,
                                1.459981699, 1.455565103, 1.451562271, 1.447940318, 1.444668631,
                                1.441718687, 1.439063886, 1.436679401, 1.434542037, 1.432630111,
                                1.430923337, 1.429402722, 1.428050477, 1.426849937, 1.425785482,
                                1.424842478, 1.424007219, 1.423266872, 1.422609439, 1.422023711,
                                1.42149924, 1.421026309, 1.420595902, 1.420199685, 1.419829987,
                                1.419479779, 1.419142664, 1.418812857, 1.41848518, 1.41815504,
                                1.417818427, 1.417471897, 1.417112559, 1.416738066, 1.416346602,
                                1.415936863, 1.415508046, 1.415059829, 1.414592355, 1.41410621,
                                1.413602403, 1.413082337, 1.412547791, 1.412000885, 1.411444054,
                                1.410880015, 1.410311734, 1.40974239, 1.409175337, 1.408614066,
                                1.408062163, 1.40752327, 1.407001038, 1.406499087, 1.406020956,
                                1.405570065, 1.405149661, 1.40476278, 1.404412198, 1.404100388,
                                1.403829476, 1.403601203, 1.403416881, 1.40327736, 1.40318299, 1.403133596, 1.403128446, 1.403166233, 1.403245055,
                                1.403362408, 1.403515176, 1.403699635};

  std::vector<std::tuple<G4double, G4double, G4double>> silGelTable;
  for (int i = 0; i < 87; ++i)
  {
    silGelTable.emplace_back(ENERGY_SilGel[i], ABSORPTION_SilGel[i], RINDEX_SilGel[i]);
  }

  // 2. Sort by energy (increasing)
  std::sort(silGelTable.begin(), silGelTable.end(),
            [](const auto &a, const auto &b)
            {
              return std::get<0>(a) < std::get<0>(b);
            });

  // 3. Remove duplicate energies (Geant4 crashes otherwise)
  auto last = std::unique(silGelTable.begin(), silGelTable.end(),
                          [](const auto &a, const auto &b)
                          {
                            return std::fabs(std::get<0>(a) - std::get<0>(b)) < 1e-9;
                          });
  silGelTable.erase(last, silGelTable.end());

  // 4. Copy back to arrays
  int N = silGelTable.size(); // CORRECT SIZE
  for (int i = 0; i < N; ++i)
  {
    ENERGY_SilGel[i] = std::get<0>(silGelTable[i]);
    ABSORPTION_SilGel[i] = std::get<1>(silGelTable[i]);
    RINDEX_SilGel[i] = std::get<2>(silGelTable[i]);
  }

  G4double RINDEX1[NUMENTRIES_water] = {1.32885, 1.32906, 1.32927, 1.32948, 1.3297, 1.32992, 1.33014,
                                        1.33037, 1.3306, 1.33084, 1.33109, 1.33134, 1.3316, 1.33186, 1.33213,
                                        1.33241, 1.3327, 1.33299, 1.33329, 1.33361, 1.33393, 1.33427, 1.33462,
                                        1.33498, 1.33536, 1.33576, 1.33617, 1.3366, 1.33705, 1.33753, 1.33803,
                                        1.33855, 1.33911, 1.3397, 1.34033, 1.341, 1.34172, 1.34248, 1.34331,
                                        1.34419, 1.34515, 1.3462, 1.34733, 1.34858, 1.34994, 1.35145, 1.35312,
                                        1.35498, 1.35707, 1.35943, 1.36211, 1.36518, 1.36872, 1.37287, 1.37776,
                                        1.38362, 1.39074, 1.39956, 1.41075, 1.42535};

  G4double ABWFF = 1.30;

  G4double ABSORPTION_water[NUMENTRIES_water] = {16.1419 * cm * ABWFF, 18.278 * cm * ABWFF, 21.0657 * cm * ABWFF, 24.8568 * cm * ABWFF,
                                                 30.3117 * cm * ABWFF, 38.8341 * cm * ABWFF, 54.0231 * cm * ABWFF, 81.2306 * cm * ABWFF,
                                                 120.909 * cm * ABWFF, 160.238 * cm * ABWFF, 193.771 * cm * ABWFF, 215.017 * cm * ABWFF,
                                                 227.747 * cm * ABWFF, 243.85 * cm * ABWFF, 294.036 * cm * ABWFF, 321.647 * cm * ABWFF,
                                                 342.81 * cm * ABWFF, 362.827 * cm * ABWFF, 378.041 * cm * ABWFF, 449.378 * cm * ABWFF,
                                                 739.434 * cm * ABWFF, 1114.23 * cm * ABWFF, 1435.56 * cm * ABWFF, 1611.06 * cm * ABWFF,
                                                 1764.18 * cm * ABWFF, 2100.95 * cm * ABWFF, 2292.9 * cm * ABWFF, 2431.33 * cm * ABWFF,
                                                 3053.6 * cm * ABWFF, 4838.23 * cm * ABWFF, 6539.65 * cm * ABWFF, 7682.63 * cm * ABWFF,
                                                 9137.28 * cm * ABWFF, 12220.9 * cm * ABWFF, 15270.7 * cm * ABWFF, 19051.5 * cm * ABWFF,
                                                 23671.3 * cm * ABWFF, 29191.1 * cm * ABWFF, 35567.9 * cm * ABWFF, 42583 * cm * ABWFF,
                                                 49779.6 * cm * ABWFF, 56465.3 * cm * ABWFF, 61830 * cm * ABWFF, 65174.6 * cm * ABWFF,
                                                 66143.7 * cm * ABWFF, 64820 * cm * ABWFF, 61635 * cm * ABWFF, 57176.2 * cm * ABWFF,
                                                 52012.1 * cm * ABWFF, 46595.7 * cm * ABWFF, 41242.1 * cm * ABWFF, 36146.3 * cm * ABWFF,
                                                 31415.4 * cm * ABWFF, 27097.8 * cm * ABWFF, 23205.7 * cm * ABWFF, 19730.3 * cm * ABWFF,
                                                 16651.6 * cm * ABWFF, 13943.6 * cm * ABWFF, 11578.1 * cm * ABWFF, 9526.13 * cm * ABWFF};

  G4double RAYFF = 0.75;
  G4double RAYLEIGH_water[NUMENTRIES_water] = {
      386929 * cm * RAYFF, 366249 * cm * RAYFF, 346398 * cm * RAYFF, 327355 * cm * RAYFF, 309097 * cm * RAYFF,
      291603 * cm * RAYFF, 274853 * cm * RAYFF, 258825 * cm * RAYFF, 243500 * cm * RAYFF, 228856 * cm * RAYFF,
      214873 * cm * RAYFF, 201533 * cm * RAYFF, 188816 * cm * RAYFF, 176702 * cm * RAYFF, 165173 * cm * RAYFF,
      154210 * cm * RAYFF, 143795 * cm * RAYFF, 133910 * cm * RAYFF, 124537 * cm * RAYFF, 115659 * cm * RAYFF,
      107258 * cm * RAYFF, 99318.2 * cm * RAYFF, 91822.2 * cm * RAYFF, 84754 * cm * RAYFF, 78097.3 * cm * RAYFF,
      71836.5 * cm * RAYFF, 65956 * cm * RAYFF, 60440.6 * cm * RAYFF, 55275.4 * cm * RAYFF, 50445.6 * cm * RAYFF,
      45937 * cm * RAYFF, 41735.2 * cm * RAYFF, 37826.6 * cm * RAYFF, 34197.6 * cm * RAYFF, 30834.9 * cm * RAYFF,
      27725.4 * cm * RAYFF, 24856.6 * cm * RAYFF, 22215.9 * cm * RAYFF, 19791.3 * cm * RAYFF, 17570.9 * cm * RAYFF,
      15543 * cm * RAYFF, 13696.6 * cm * RAYFF, 12020.5 * cm * RAYFF, 10504.1 * cm * RAYFF, 9137.15 * cm * RAYFF,
      7909.45 * cm * RAYFF, 6811.3 * cm * RAYFF, 5833.25 * cm * RAYFF, 4966.2 * cm * RAYFF, 4201.36 * cm * RAYFF,
      3530.28 * cm * RAYFF, 2944.84 * cm * RAYFF, 2437.28 * cm * RAYFF, 2000.18 * cm * RAYFF, 1626.5 * cm * RAYFF,
      1309.55 * cm * RAYFF, 1043.03 * cm * RAYFF, 821.016 * cm * RAYFF, 637.97 * cm * RAYFF, 488.754 * cm * RAYFF};

  G4double ENERGY_skAcrylic[306] =
      {2.066 * eV, 2.070 * eV, 2.073 * eV, 2.077 * eV, 2.080 * eV, 2.084 * eV, 2.087 * eV,
       2.091 * eV, 2.094 * eV, 2.098 * eV, 2.101 * eV, 2.105 * eV, 2.109 * eV, 2.112 * eV,
       2.116 * eV, 2.119 * eV, 2.123 * eV, 2.127 * eV, 2.130 * eV, 2.134 * eV, 2.138 * eV,
       2.141 * eV, 2.145 * eV, 2.149 * eV, 2.153 * eV, 2.156 * eV, 2.160 * eV, 2.164 * eV,
       2.168 * eV, 2.171 * eV, 2.175 * eV, 2.179 * eV, 2.183 * eV, 2.187 * eV, 2.191 * eV,
       2.194 * eV, 2.198 * eV, 2.202 * eV, 2.206 * eV, 2.210 * eV, 2.214 * eV, 2.218 * eV,
       2.222 * eV, 2.226 * eV, 2.230 * eV, 2.234 * eV, 2.238 * eV, 2.242 * eV, 2.246 * eV,
       2.250 * eV, 2.254 * eV, 2.258 * eV, 2.262 * eV, 2.267 * eV, 2.271 * eV, 2.275 * eV,
       2.279 * eV, 2.283 * eV, 2.288 * eV, 2.292 * eV, 2.296 * eV, 2.300 * eV, 2.305 * eV,
       2.309 * eV, 2.313 * eV, 2.317 * eV, 2.322 * eV, 2.326 * eV, 2.331 * eV, 2.335 * eV,
       2.339 * eV, 2.344 * eV, 2.348 * eV, 2.353 * eV, 2.357 * eV, 2.362 * eV, 2.366 * eV,
       2.371 * eV, 2.375 * eV, 2.380 * eV, 2.384 * eV, 2.389 * eV, 2.394 * eV, 2.398 * eV,
       2.403 * eV, 2.407 * eV, 2.412 * eV, 2.417 * eV, 2.422 * eV, 2.426 * eV, 2.431 * eV,
       2.436 * eV, 2.441 * eV, 2.445 * eV, 2.450 * eV, 2.455 * eV, 2.460 * eV, 2.465 * eV,
       2.470 * eV, 2.475 * eV, 2.480 * eV, 2.485 * eV, 2.490 * eV, 2.495 * eV, 2.500 * eV,
       2.505 * eV, 2.510 * eV, 2.515 * eV, 2.520 * eV, 2.525 * eV, 2.530 * eV, 2.535 * eV,
       2.541 * eV, 2.546 * eV, 2.551 * eV, 2.556 * eV, 2.562 * eV, 2.567 * eV, 2.572 * eV,
       2.578 * eV, 2.583 * eV, 2.588 * eV, 2.594 * eV, 2.599 * eV, 2.605 * eV, 2.610 * eV,
       2.616 * eV, 2.621 * eV, 2.627 * eV, 2.632 * eV, 2.638 * eV, 2.644 * eV, 2.649 * eV,
       2.655 * eV, 2.661 * eV, 2.666 * eV, 2.672 * eV, 2.678 * eV, 2.684 * eV, 2.689 * eV,
       2.695 * eV, 2.701 * eV, 2.707 * eV, 2.713 * eV, 2.719 * eV, 2.725 * eV, 2.731 * eV,
       2.737 * eV, 2.743 * eV, 2.749 * eV, 2.755 * eV, 2.761 * eV, 2.768 * eV, 2.774 * eV,
       2.780 * eV, 2.786 * eV, 2.792 * eV, 2.799 * eV, 2.805 * eV, 2.811 * eV, 2.818 * eV,
       2.824 * eV, 2.831 * eV, 2.837 * eV, 2.844 * eV, 2.850 * eV, 2.857 * eV, 2.863 * eV,
       2.870 * eV, 2.877 * eV, 2.883 * eV, 2.890 * eV, 2.897 * eV, 2.904 * eV, 2.910 * eV,
       2.917 * eV, 2.924 * eV, 2.931 * eV, 2.938 * eV, 2.945 * eV, 2.952 * eV, 2.959 * eV,
       2.966 * eV, 2.973 * eV, 2.980 * eV, 2.988 * eV, 2.995 * eV, 3.002 * eV, 3.009 * eV,
       3.017 * eV, 3.024 * eV, 3.031 * eV, 3.039 * eV, 3.046 * eV, 3.054 * eV, 3.061 * eV,
       3.069 * eV, 3.077 * eV, 3.084 * eV, 3.092 * eV, 3.100 * eV, 3.107 * eV, 3.115 * eV,
       3.123 * eV, 3.131 * eV, 3.139 * eV, 3.147 * eV, 3.155 * eV, 3.163 * eV, 3.171 * eV,
       3.179 * eV, 3.187 * eV, 3.195 * eV, 3.204 * eV, 3.212 * eV, 3.220 * eV, 3.229 * eV,
       3.237 * eV, 3.246 * eV, 3.254 * eV, 3.263 * eV, 3.271 * eV, 3.280 * eV, 3.289 * eV,
       3.297 * eV, 3.306 * eV, 3.315 * eV, 3.324 * eV, 3.333 * eV, 3.342 * eV, 3.351 * eV,
       3.360 * eV, 3.369 * eV, 3.378 * eV, 3.388 * eV, 3.397 * eV, 3.406 * eV, 3.416 * eV,
       3.425 * eV, 3.434 * eV, 3.444 * eV, 3.454 * eV, 3.463 * eV, 3.473 * eV, 3.483 * eV,
       3.493 * eV, 3.502 * eV, 3.512 * eV, 3.522 * eV, 3.532 * eV, 3.542 * eV, 3.553 * eV,
       3.563 * eV, 3.573 * eV, 3.583 * eV, 3.594 * eV, 3.604 * eV, 3.615 * eV, 3.625 * eV,
       3.636 * eV, 3.647 * eV, 3.657 * eV, 3.668 * eV, 3.679 * eV, 3.690 * eV, 3.701 * eV,
       3.712 * eV, 3.723 * eV, 3.734 * eV, 3.746 * eV, 3.757 * eV, 3.769 * eV, 3.780 * eV,
       3.792 * eV, 3.803 * eV, 3.815 * eV, 3.827 * eV, 3.839 * eV, 3.850 * eV, 3.862 * eV,
       3.875 * eV, 3.887 * eV, 3.899 * eV, 3.911 * eV, 3.924 * eV, 3.936 * eV, 3.949 * eV,
       3.961 * eV, 3.974 * eV, 3.987 * eV, 3.999 * eV, 4.012 * eV, 4.025 * eV, 4.039 * eV,
       4.052 * eV, 4.065 * eV, 4.078 * eV, 4.092 * eV, 4.105 * eV, 4.119 * eV, 4.133 * eV,
       4.147 * eV, 4.161 * eV, 4.175 * eV, 4.189 * eV, 4.203 * eV};

  G4double RINDEX_skAcrylic[306] =
      {1.4901, 1.4901, 1.4902, 1.4902, 1.4902, 1.4903, 1.4903, 1.4904, 1.4904, 1.4904,
       1.4905, 1.4905, 1.4906, 1.4906, 1.4907, 1.4907, 1.4908, 1.4908, 1.4908, 1.4909,
       1.4909, 1.4910, 1.4910, 1.4911, 1.4911, 1.4912, 1.4912, 1.4913, 1.4913, 1.4913,
       1.4914, 1.4914, 1.4915, 1.4915, 1.4916, 1.4916, 1.4917, 1.4917, 1.4918, 1.4918,
       1.4919, 1.4919, 1.4920, 1.4920, 1.4921, 1.4921, 1.4922, 1.4922, 1.4923, 1.4923,
       1.4924, 1.4924, 1.4925, 1.4926, 1.4926, 1.4927, 1.4927, 1.4928, 1.4928, 1.4929,
       1.4929, 1.4930, 1.4931, 1.4931, 1.4932, 1.4932, 1.4933, 1.4933, 1.4934, 1.4935,
       1.4935, 1.4936, 1.4936, 1.4937, 1.4938, 1.4938, 1.4939, 1.4940, 1.4940, 1.4941,
       1.4941, 1.4942, 1.4943, 1.4943, 1.4944, 1.4945, 1.4945, 1.4946, 1.4947, 1.4947,
       1.4948, 1.4949, 1.4949, 1.4950, 1.4951, 1.4951, 1.4952, 1.4953, 1.4954, 1.4954,
       1.4955, 1.4956, 1.4956, 1.4957, 1.4958, 1.4959, 1.4959, 1.4960, 1.4961, 1.4962,
       1.4962, 1.4963, 1.4964, 1.4965, 1.4965, 1.4966, 1.4967, 1.4968, 1.4969, 1.4970,
       1.4970, 1.4971, 1.4972, 1.4973, 1.4974, 1.4974, 1.4975, 1.4976, 1.4977, 1.4978,
       1.4979, 1.4980, 1.4980, 1.4981, 1.4982, 1.4983, 1.4984, 1.4985, 1.4986, 1.4987,
       1.4988, 1.4989, 1.4990, 1.4991, 1.4992, 1.4993, 1.4994, 1.4995, 1.4996, 1.4997,
       1.4997, 1.4998, 1.5000, 1.5001, 1.5002, 1.5003, 1.5004, 1.5005, 1.5006, 1.5007,
       1.5008, 1.5009, 1.5010, 1.5011, 1.5012, 1.5013, 1.5015, 1.5016, 1.5017, 1.5018,
       1.5019, 1.5020, 1.5022, 1.5023, 1.5024, 1.5025, 1.5026, 1.5028, 1.5029, 1.5030,
       1.5031, 1.5033, 1.5034, 1.5035, 1.5036, 1.5038, 1.5039, 1.5040, 1.5042, 1.5043,
       1.5044, 1.5046, 1.5047, 1.5049, 1.5050, 1.5051, 1.5053, 1.5054, 1.5056, 1.5057,
       1.5059, 1.5060, 1.5062, 1.5063, 1.5065, 1.5066, 1.5068, 1.5069, 1.5071, 1.5073,
       1.5074, 1.5076, 1.5077, 1.5079, 1.5081, 1.5082, 1.5084, 1.5086, 1.5088, 1.5089,
       1.5091, 1.5093, 1.5095, 1.5096, 1.5098, 1.5100, 1.5102, 1.5104, 1.5106, 1.5108,
       1.5109, 1.5111, 1.5113, 1.5115, 1.5117, 1.5119, 1.5121, 1.5123, 1.5125, 1.5127,
       1.5130, 1.5132, 1.5134, 1.5136, 1.5138, 1.5141, 1.5143, 1.5145, 1.5147, 1.5149,
       1.5152, 1.5154, 1.5157, 1.5159, 1.5161, 1.5164, 1.5166, 1.5169, 1.5171, 1.5174,
       1.5176, 1.5179, 1.5181, 1.5184, 1.5187, 1.5190, 1.5192, 1.5195, 1.5198, 1.5201,
       1.5203, 1.5206, 1.5209, 1.5212, 1.5215, 1.5218, 1.5221, 1.5224, 1.5227, 1.5230,
       1.5234, 1.5237, 1.5240, 1.5243, 1.5247, 1.5250, 1.5254, 1.5257, 1.5260, 1.5264,
       1.5267, 1.5271, 1.5275, 1.5279, 1.5282, 1.5286, 1.5290, 1.5294, 1.5297, 1.5302,
       1.5306, 1.5310, 1.5314, 1.5318, 1.5322, 1.5326};

  G4double ABSORPTION_skAcrylic[306] =
      {3472.22 * m, 3355.70 * m, 3300.33 * m, 3412.97 * m, 3496.50 * m, 3649.64 * m, 3861.00 * m, 3921.57 * m, 3846.15 * m, 3703.70 * m,
       3448.28 * m, 3048.78 * m, 2932.55 * m, 2770.08 * m, 2808.99 * m, 2906.98 * m, 2898.55 * m, 2906.98 * m, 3003.00 * m, 3164.56 * m,
       3184.71 * m, 3246.75 * m, 3086.42 * m, 2994.01 * m, 2994.01 * m, 2967.36 * m, 3030.30 * m, 3039.51 * m, 3105.59 * m, 2958.58 * m,
       3012.05 * m, 3039.51 * m, 2967.36 * m, 2932.55 * m, 2923.98 * m, 2890.17 * m, 2857.14 * m, 2865.33 * m, 2915.45 * m, 2898.55 * m,
       2881.84 * m, 2754.82 * m, 2702.70 * m, 2754.82 * m, 2890.17 * m, 2840.91 * m, 2659.57 * m, 2604.17 * m, 2624.67 * m, 2857.14 * m,
       2849.00 * m, 2849.00 * m, 2610.97 * m, 2597.40 * m, 2688.17 * m, 2590.67 * m, 2439.02 * m, 2380.95 * m, 2398.08 * m, 2320.19 * m,
       2369.67 * m, 2293.58 * m, 2267.57 * m, 2262.44 * m, 2320.19 * m, 2309.47 * m, 2237.14 * m, 2247.19 * m, 2257.34 * m, 2336.45 * m,
       2347.42 * m, 2421.31 * m, 2439.02 * m, 2392.34 * m, 2392.34 * m, 2500.00 * m, 2500.00 * m, 2512.56 * m, 2481.39 * m, 2444.99 * m,
       2155.17 * m, 2053.39 * m, 2008.03 * m, 2040.82 * m, 1988.07 * m, 2004.01 * m, 1976.28 * m, 1834.86 * m, 1798.56 * m, 1869.16 * m,
       1808.32 * m, 1754.39 * m, 1805.05 * m, 1841.62 * m, 1872.66 * m, 1964.64 * m, 1976.28 * m, 1893.94 * m, 1890.36 * m, 1908.40 * m,
       1897.53 * m, 1855.29 * m, 1886.79 * m, 1904.76 * m, 1841.62 * m, 1855.29 * m, 1834.86 * m, 1805.05 * m, 1733.10 * m, 1779.36 * m,
       1754.39 * m, 1788.91 * m, 1831.50 * m, 1834.86 * m, 1869.16 * m, 1912.05 * m, 2004.01 * m, 1968.50 * m, 1915.71 * m, 1818.18 * m,
       1795.33 * m, 1766.78 * m, 1808.32 * m, 1766.78 * m, 1754.39 * m, 1751.31 * m, 1776.20 * m, 1838.24 * m, 1872.66 * m, 1919.39 * m,
       1901.14 * m, 1988.07 * m, 1926.78 * m, 1893.94 * m, 1901.14 * m, 1851.85 * m, 1818.18 * m, 1769.91 * m, 1727.12 * m, 1724.14 * m,
       1709.40 * m, 1689.19 * m, 1692.05 * m, 1680.67 * m, 1697.79 * m, 1727.12 * m, 1763.67 * m, 1697.79 * m, 1692.05 * m, 1666.67 * m,
       1600.00 * m, 1636.66 * m, 1615.51 * m, 1610.31 * m, 1538.46 * m, 1501.50 * m, 1519.76 * m, 1492.54 * m, 1515.15 * m, 1457.73 * m,
       1479.29 * m, 1461.99 * m, 1470.59 * m, 1470.59 * m, 1453.49 * m, 1512.86 * m, 1522.07 * m, 1589.83 * m, 1577.29 * m, 1597.44 * m,
       1600.00 * m, 1672.24 * m, 1663.89 * m, 1647.45 * m, 1661.13 * m, 1607.72 * m, 1543.21 * m, 1488.10 * m, 1459.85 * m, 1373.63 * m,
       1326.26 * m, 1373.63 * m, 1410.44 * m, 1406.47 * m, 1396.65 * m, 1308.90 * m, 1328.02 * m, 1386.96 * m, 1347.71 * m, 1297.02 * m,
       1262.63 * m, 1234.57 * m, 1283.70 * m, 1319.26 * m, 1310.62 * m, 1319.26 * m, 1402.52 * m, 1459.85 * m, 1464.13 * m, 1472.75 * m,
       1420.45 * m, 1432.66 * m, 1400.56 * m, 1362.40 * m, 1358.70 * m, 1358.70 * m, 1379.31 * m, 1360.54 * m, 1362.40 * m, 1355.01 * m,
       1342.28 * m, 1367.99 * m, 1293.66 * m, 1221.00 * m, 1197.60 * m, 1165.50 * m, 1119.82 * m, 1107.42 * m, 1108.65 * m, 1040.58 * m,
       1033.06 * m, 1000.00 * m, 952.38 * m, 925.93 * m, 909.09 * m, 869.57 * m, 800.00 * m, 781.25 * m, 769.23 * m, 729.93 * m, 735.29 * m,
       729.93 * m, 709.22 * m, 684.93 * m, 657.89 * m, 645.16 * m, 625.00 * m, 581.40 * m, 568.18 * m, 558.66 * m, 543.48 * m, 540.54 * m,
       523.56 * m, 505.05 * m, 490.20 * m, 485.44 * m, 462.96 * m, 450.45 * m, 442.48 * m, 429.18 * m, 420.17 * m, 411.52 * m, 411.52 * m,
       393.70 * m, 393.70 * m, 386.10 * m, 381.68 * m, 355.87 * m, 343.64 * m, 331.13 * m, 343.64 * m, 338.98 * m, 335.57 * m, 326.80 * m,
       313.48 * m, 303.95 * m, 294.12 * m, 292.40 * m, 280.90 * m, 273.97 * m, 270.27 * m, 265.25 * m, 261.78 * m, 256.41 * m, 250.00 * m,
       242.13 * m, 237.53 * m, 233.10 * m, 228.83 * m, 226.76 * m, 223.21 * m, 219.30 * m, 215.05 * m, 211.86 * m, 208.77 * m, 206.61 * m,
       201.21 * m, 196.85 * m, 193.42 * m, 189.75 * m, 186.22 * m, 182.15 * m, 175.75 * m, 168.92 * m, 162.60 * m, 154.80 * m, 144.51 * m,
       133.51 * m, 120.19 * m, 105.93 * m, 90.91 * m, 75.76 * m, 61.35 * m, 48.54 * m, 37.88 * m, 29.41 * m};

  G4double Matrix_ref[NUMENTRIES_water] =
      {.0, .0, .0, .0, .0, .0,
       .0, .0, .0, .0, .0, .0,
       .0, .0, .0, .0, .0, .0,
       .0, .0, .0, .0, .0, .0,
       .0, .0, .0, .0, .0, .0,
       .0, .0, .0, .0, .0, .0,
       .0, .0, .0, .0, .0, .0,
       .0, .0, .0, .0, .0, .0,
       .0, .0, .0, .0, .0, .0,
       .0, .0, .0, .0, .0, .0};

  ////////////AIR//////////////////////
  a = 14.01 * g / mole;
  G4Element *elN = new G4Element("Nitrogen", "N", 7, a);

  a = 16.00 * g / mole;
  G4Element *elO = new G4Element("Oxygen", "O", 8, a);

  density = 1.290 * mg / cm3;
  G4Material *Air = new G4Material("Air", density, 2);
  Air->AddElement(elN, 70. * perCent);
  Air->AddElement(elO, 30. * perCent);

  G4MaterialPropertiesTable *myMPT2 = new G4MaterialPropertiesTable();
  myMPT2->AddProperty("RINDEX", ENERGY_water, RINDEX_air, NUMENTRIES_water);

  Air->SetMaterialPropertiesTable(myMPT2);

  ///////////////////////Silica Gel///////////////////////
  a = 28.09 * g / mole;
  G4Element *elSi = new G4Element("Silicon", "Si", 14., a);

  density = 0.98 * g / cm3;
  G4Material *SilGel = new G4Material("SilGel", density, 1);
  SilGel->AddElement(elSi, 1);

  G4MaterialPropertiesTable *SilGelPropTable = new G4MaterialPropertiesTable();
  SilGelPropTable->AddProperty("RINDEX", ENERGY_SilGel, RINDEX_SilGel, 87);
  SilGelPropTable->AddProperty("ABSLENGTH", ENERGY_SilGel, ABSORPTION_SilGel, 87);
  SilGel->SetMaterialPropertiesTable(SilGelPropTable);

  ////////////////////PLA///////////////////////

  a = 1.01 * g / mole;
  G4Element *elH = new G4Element("Hydrogen", "H", 1, a);
  G4Element *elC = new G4Element("Carbon", "C", 6, 12.01 * g / mole);

  density = 0.95 * g / cm3;
  G4Material *Plastic = new G4Material("Plastic", density, 2);
  Plastic->AddElement(elC, 1);
  Plastic->AddElement(elH, 2);

  G4MaterialPropertiesTable *myMPT3 = new G4MaterialPropertiesTable();
  myMPT3->AddProperty("ABSLENGTH", ENERGY_water, Matrix_ref, NUMENTRIES_water);
  Plastic->SetMaterialPropertiesTable(myMPT3);

  G4OpticalSurface *PlasticSkinSurface = new G4OpticalSurface("PlasticSurface");
  PlasticSkinSurface->SetType(dielectric_dielectric);
  PlasticSkinSurface->SetModel(unified);
  PlasticSkinSurface->SetFinish(polished);

  //////////////Water///////////////////
  density = 1.00 * g / cm3;
  G4Material *Water = new G4Material("Water", density, 2);
  Water->AddElement(elH, 2);
  Water->AddElement(elO, 1);

  G4MaterialPropertiesTable *myMPT1 = new G4MaterialPropertiesTable();
  myMPT1->AddProperty("RINDEX", ENERGY_water, RINDEX1, NUMENTRIES_water);
  myMPT1->AddProperty("ABSLENGTH", ENERGY_water, ABSORPTION_water, NUMENTRIES_water);
  myMPT1->AddProperty("RAYLEIGH", ENERGY_water, RAYLEIGH_water, NUMENTRIES_water);
  Water->SetMaterialPropertiesTable(myMPT1);

  ////////////////////////////Acrylic//////////////////////
  G4Material *Acrylic = nist->FindOrBuildMaterial("G4_PLEXIGLASS");

  G4MaterialPropertiesTable *AcrPropTable = new G4MaterialPropertiesTable();
  AcrPropTable->AddProperty("RINDEX", ENERGY_skAcrylic, RINDEX_skAcrylic, 306);
  AcrPropTable->AddProperty("ABSLENGTH", ENERGY_skAcrylic, ABSORPTION_skAcrylic, 306);
  Acrylic->SetMaterialPropertiesTable(AcrPropTable);

  ////////////////////Making holes in the Matrix////////////////

  WCSimTuningParameters *tuningpars = new WCSimTuningParameters();
  enum DetConfiguration
  {
    wfm = 1,
    fwm = 2
  };
  G4int WCSimConfiguration = fwm;
  inSituPMTConstruction inSituPMT(WCSimConfiguration, tuningpars);

  G4int numPMTs = 19;

  G4double vessel_cylinder_height = 77.785 * 2 * mm;
  G4double vessel_inner_radius = 250.93 * mm;
  G4double vessel_outer_radius = 254. * mm;

  G4double domeInnerRadius = 332. * mm;
  G4double domeOuterRadius = 347. * mm;
  G4double domeCut = 235 * mm;
  G4double dome_height = domeOuterRadius - domeCut;

  G4double mPMT_zRange_outer[2] = {0, // start from zero for easier placement in ID
                                   vessel_cylinder_height};
  G4double mPMT_RRange_outer[2] = {vessel_outer_radius,
                                   vessel_outer_radius};
  G4double mPMT_rRange_outer[2] = {0., 0.};

  // Although G4Tubs is more natural, Polycone is used to be in control of z position
  // and because z = 0 is position of bottom of cylinder, vs center in G4Tubs
  G4Polycone *solidMultiPMT_vessel =
      new G4Polycone("WCmPMT_vessel",
                     0.0 * deg,
                     360.0 * deg,
                     2,
                     mPMT_zRange_outer,
                     mPMT_rRange_outer,  // R Inner
                     mPMT_RRange_outer); // R Outer

  G4Sphere *mPMT_top_sphere =
      new G4Sphere("WCmPMT_tsphere",
                   0,
                   domeOuterRadius,
                   0.0 * deg, 360.0 * deg,
                   0.0 * deg, 90.0 * deg);

  G4Box *solidCutOffTubs =
      new G4Box("cutOffTubs",
                domeOuterRadius + 1. * cm, // x/2
                domeOuterRadius + 1. * cm, // y/2
                domeCut);                  // z/2

  G4SubtractionSolid *mPMT_top_cap_vessel = new G4SubtractionSolid("WCmPMT_tsphere_vessel",
                                                                   mPMT_top_sphere,
                                                                   solidCutOffTubs);

  G4UnionSolid *solidMultiPMT =
      new G4UnionSolid("WCmPMT", solidMultiPMT_vessel, mPMT_top_cap_vessel, 0, G4ThreeVector(0, 0, vessel_cylinder_height - domeCut));

  G4Box *solidInSitumPMT = new G4Box("solidInSitumPMT", 3.5 * m, 3.5 * m, 3.5 * m);
  G4LogicalVolume *logicWCMultiPMT_outer =
      new G4LogicalVolume(solidInSitumPMT,
                          Air,
                          // Water, 
                          "WCMultiPMT_outer",
                          0, 0, 0);

  G4VPhysicalVolume *outerWorld = new G4PVPlacement(0,
                                                   G4ThreeVector(0., 0., 0.),
                                                   logicWCMultiPMT_outer,
                                                   "outerWorld",
                                                   0,
                                                   false,
                                                   0,
                                                   true);

  G4Box *WaterLevel = new G4Box("WaterLevel", 3.49 * m, 3.49 * m, 3.49 * m);
  G4LogicalVolume *logicWCMultiPMT =
      new G4LogicalVolume(WaterLevel,
                          //Air,
                          Water, 
                          "WCMultiPMT",
                          0, 0, 0);       
                          
  G4VPhysicalVolume *physWorld = new G4PVPlacement(0,
                                                   G4ThreeVector(0., 0., 0.),
                                                   logicWCMultiPMT,
                                                   "physWorld",
                                                   logicWCMultiPMT_outer,
                                                   false,
                                                   0,
                                                   true);

  G4VisAttributes *worldAttributes = new G4VisAttributes();
  worldAttributes->SetColor(0.5, 0.5, 0.5, 0.5);
  worldAttributes->SetVisibility(true);
  worldAttributes->SetForceSolid(true);
  logicWCMultiPMT->SetVisAttributes(worldAttributes);

  G4Sphere *domeSphere =
      new G4Sphere("DomeSphere",
                   domeInnerRadius,
                   domeOuterRadius,
                   0.0 * deg, 360.0 * deg,
                   0.0, 90. * deg);

  G4VSolid *domeSolid = new G4SubtractionSolid("domeSolid",
                                               domeSphere,
                                               solidCutOffTubs);

  G4LogicalVolume *domeLogic = new G4LogicalVolume(domeSolid,
                                                   Acrylic,
                                                   "logicDome");

  new G4PVPlacement(0,
                    G4ThreeVector(0., 0., vessel_cylinder_height - domeCut),
                    domeLogic,
                    "physInsituDome",
                    logicWCMultiPMT,
                    false,
                    0,
                    true);

  G4VisAttributes *domeAttributes = new G4VisAttributes();
  domeAttributes->SetColor(1.0, 1.0, 1.0, 0.5);
  domeAttributes->SetVisibility(true);
  domeAttributes->SetForceSolid(true);
  domeLogic->SetVisAttributes(domeAttributes);

  G4Tubs *cylinderSolid = new G4Tubs("cylinderSolid",
                                     vessel_inner_radius,
                                     vessel_outer_radius,
                                     vessel_cylinder_height / 2,
                                     0. * deg, 360. * deg);

  G4LogicalVolume *logicCylinder = new G4LogicalVolume(cylinderSolid,
                                                       Plastic,
                                                       "logicCylinder");

  new G4PVPlacement(0,
                    G4ThreeVector(0., 0., vessel_cylinder_height / 2),
                    logicCylinder,
                    "physInsituCylinder",
                    logicWCMultiPMT,
                    false,
                    0,
                    true);

  new G4LogicalSkinSurface("cylinderSkin", logicCylinder, PlasticSkinSurface);

  G4VisAttributes *cylinderAttributes = new G4VisAttributes();
  cylinderAttributes->SetColor(0.2, 0.2, 0.2, 1.0);
  cylinderAttributes->SetVisibility(true);
  cylinderAttributes->SetForceSolid(true);
  logicCylinder->SetVisAttributes(cylinderAttributes);

  /////////////////////////////////////////////////////////////////////////
  /// 3) This is the area between the outer shell and the inner shell.  ///
  ///    In this space, the PMTs will live.                             ///
  /////////////////////////////////////////////////////////////////////////

  // Matrix, where PMTs are placed as daughters
  G4double matrixRmin = 263.003 * mm;
  G4double matrixRmax = 325.603 * mm;
  G4double matrixUpshift = 5.7 * mm;
  G4Sphere *solidMatrix =
      new G4Sphere("solidMatrix",
                   matrixRmin,
                   matrixRmax,
                   0.0 * deg, 360.0 * deg,
                   0.0, 43.4549 * deg);

  G4LogicalVolume *logicMatrix = new G4LogicalVolume(solidMatrix,
                                                     Plastic,
                                                     "logicMatrix");

  new G4LogicalSkinSurface("matrixSkinSurface", logicMatrix, PlasticSkinSurface);

  G4double gelRmin = 323.5 * mm;
  G4Sphere *solidGelSphere = new G4Sphere("solidGelSphere",
                                          gelRmin,
                                          domeInnerRadius,
                                          0.0 * deg, 360.0 * deg,
                                          0.0, 90.0 * deg);

  G4Box *solidGelCutOut = new G4Box("BoxGelCutOut",
                                    domeInnerRadius + 1. * cm,
                                    domeInnerRadius + 1. * cm,
                                    241.04 * mm);

  G4VSolid *solidGelSphereCut =
      new G4SubtractionSolid("solidGelSphereCut",
                             solidGelSphere,
                             solidGelCutOut);

  G4VSolid *solidGel = new G4SubtractionSolid("solidGel",
                                              solidGelSphereCut,
                                              solidMatrix,
                                              0,
                                              G4ThreeVector(0., 0., gelRmin - matrixRmax + matrixUpshift));

  G4LogicalVolume *logicGel = new G4LogicalVolume(solidGel,
                                                  SilGel,
                                                  "logicGel");

  new G4PVPlacement(0,
                    G4ThreeVector(0., 0., vessel_cylinder_height - domeCut),
                    logicGel,
                    "gel",
                    logicWCMultiPMT,
                    false,
                    0,
                    true);

  G4VisAttributes *gelAttributes = new G4VisAttributes();
  gelAttributes->SetColor(1.0, 1.0, 1.0, 0.5);
  gelAttributes->SetVisibility(true);
  gelAttributes->SetForceSolid(true);
  logicGel->SetVisAttributes(gelAttributes);

  // Seal the bottom with blacksheet
  G4Tubs *solidBS = new G4Tubs("solidBS",
                               0 * mm,
                               vessel_inner_radius,
                               1. * mm,
                               0. * deg, 360. * deg);

  G4LogicalVolume *logicBS = new G4LogicalVolume(solidBS,
                                                 Plastic,
                                                 "logicBS");

  new G4PVPlacement(0,
                    G4ThreeVector(0., 0., 1. * mm),
                    logicBS,
                    "physInsituMPMTBS",
                    logicWCMultiPMT,
                    false,
                    0,
                    true);

  new G4LogicalSkinSurface("InsitumPMTBSSkinSurface", logicBS, PlasticSkinSurface);

  G4VisAttributes *VisAttRed = new G4VisAttributes(G4Colour(1.0, 0., 0.));
  VisAttRed->SetForceSolid(true);
  logicBS->SetVisAttributes(VisAttRed);

  G4LogicalVolume *inSituPMTLogic = inSituPMT.ConstructioninSituPMT();
  // G4LogicalVolume* logicWCPMT = ConstructInSituPMT(PMTName, CollectionName,detectorElement);

  G4double thetaArray[19] = {0. * deg, 18.1345 * deg, 18.1345 * deg, 18.1345 * deg, 18.1345 * deg, 18.1345 * deg, 18.1345 * deg,
                             35.3149 * deg, 35.3149 * deg, 35.3149 * deg, 35.3149 * deg, 35.3149 * deg, 35.3149 * deg,
                             35.3149 * deg, 35.3149 * deg, 35.3149 * deg, 35.3149 * deg, 35.3149 * deg, 35.3149 * deg};
  // G4double phiArray[19] = {0. * deg, 0. * deg, 60. * deg, 120. * deg, 180. * deg, -120. * deg, -60 * deg, 0. * deg, 30. * deg, 60. * deg, 90. * deg, 120. * deg, 150. * deg, 180. * deg, -150. * deg, -120. * deg, -90. * deg, -60. * deg, -30. * deg};
  G4double phiArray[] = {0. * deg, 0. * deg, -60. * deg, -120. * deg, 180. * deg, 120. * deg, 60 * deg, 0. * deg, -30. * deg, -60. * deg, -90. * deg, -120. * deg, -150. * deg, 180. * deg, 150. * deg, 120. * deg, 90. * deg, 60. * deg, 30. * deg};

  // Place the PMTs inside matrix
  for (int i = 0; i < numPMTs; i++)
  {

    std::stringstream ss;
    ss << "PMTCopy" << i;
    std::string pmtName = ss.str();

    // G4ThreeVector PMTPosition = {xArray[i], yArray[i], zArray[i]};
    G4ThreeVector PMTPosition = {0.0, 0.0, 0.0};

    G4RotationMatrix *PMTRotation = new G4RotationMatrix;
    PMTRotation->rotateZ(-phiArray[i]);
    PMTRotation->rotateY(-thetaArray[i]);

    // Create and place the PMT copy
    new G4PVPlacement(PMTRotation,
                      PMTPosition,
                      inSituPMTLogic,
                      pmtName,
                      logicMatrix,
                      false,
                      i,
                      true);
  }

  // Finally place the matrix
  new G4PVPlacement(0,
                    G4ThreeVector(0., 0., vessel_cylinder_height - domeCut + gelRmin - matrixRmax + matrixUpshift),
                    logicMatrix,
                    "Matrix",
                    logicWCMultiPMT,
                    false,
                    0,
                    true);

  /* Set all visualization here for better overview. */
  // Gray wireframe visual style
  G4VisAttributes *VisAttGrey = new G4VisAttributes(G4Colour(0.2, 0.2, 0.2));
  // VisAttGrey->SetForceWireframe(true);
  VisAttGrey->SetForceSolid(true);
  G4VisAttributes *WCPMTVisAtt5 = new G4VisAttributes(G4Colour(.0, 1., 1.));
  WCPMTVisAtt5->SetForceSolid(true);
  G4VisAttributes *VisAttYellow = new G4VisAttributes(G4Colour(1.0, 1., 0.));
  VisAttYellow->SetForceSolid(true);

  G4RunManager::GetRunManager()->DefineWorldVolume(outerWorld);

  return outerWorld;
}

# Where the numbers of the made-up stars come from

The stars, companions and planets that `universe.json` has computed around the camera are not
real ones: they follow the statistics of real ones. The tables of the map were checked against
these sources on 2026-11-15. "Not verified" means no page stating the number was opened.

## Stars (`starField.brightness`, `starField.whiteDwarfShare`)

| What | Value used | Source |
|---|---|---|
| Share of main-sequence stars by class near the Sun | M 76%, K 12%, G 7.6%, F 3.0%, A 0.61%, B 0.12%, O 0.00003% | Ledrew 2001, as tabulated in https://en.wikipedia.org/wiki/Stellar_classification (the paper itself was not opened) |
| Absolute V magnitude where the classes meet | B0V -3.90, A0V 0.99, F0V 2.57, G0V 4.48, K0V 5.78, M0V 8.80, M9.5V 19.75 | Pecaut & Mamajek, https://www.pas.rochester.edu/~emamajek/EEM_dwarf_UBVIJHK_colors_Teff.txt (version 2022.04.16) |
| Mass at those points, in solar masses | 17.7, 2.18, 1.61, 1.06, 0.88, 0.57, 0.078 | same table |
| Radius and temperature along the main sequence | 27 rows of the table, O5V (11.45 solar radii, 41400 K) to M9.5V (0.101 solar radii, 2350 K); the Sun's class G2V: 1.012 radii, 5770 K | same table, read on 2026-11-15; the rows are in `viewer/src/stars.rs` |
| Giants and white dwarfs | a giant's radius follows from its brightness and a temperature of 4000 to 5200 K; a white dwarf is 0.012 solar radii, about the size of the Earth | not from a table: the usual textbook values |
| White dwarfs | about 1e10 in a galaxy of 1e11 to 4e11 stars: 2.5% of 4e11 | https://en.wikipedia.org/wiki/Milky_Way (primary references not identified) |
| Stars in all | 1e11 to 4e11; the map takes the upper end | same page |

## Companions

| What | Value | Source |
|---|---|---|
| Sun-like stars: single / double / triple / more | 56 / 33 / 8 / 3 % | Raghavan et al. 2010, https://arxiv.org/abs/1007.0414 |
| Sun-like stars: periods | log-normal, mean log P(days) = 5.03, sigma 2.28 (about 45 AU) | same |
| M dwarfs: single / double / triple / more | 76 / 20 / 3 / 0.3 % | Winters et al. 2019, https://arxiv.org/abs/1901.06364 |
| M dwarfs: separations | log-normal, peak 20 AU, sigma(log a) = 1.16 (Duchene & Kraus 2013 give 5.3 AU) | same; https://arxiv.org/abs/1303.3028 |
| Multiple / triple-or-more by mass | early M 30 / 6.3 %, FGK 46 / 12 %, A 68 / 25 %, B (3-5 Msun) 81 / 36 %, B (8-17) 93 / 57 %, O 96 / 68 % | Offner et al. 2023, table 1, https://arxiv.org/abs/2203.10066 |
| Mass ratios | roughly flat from 0.2 to 0.95 for sun-like stars, 0.5 to 1 for M dwarfs | Raghavan 2010; Winters 2019 |
| Triples that hold together | a_out / a_in > 2.8 / (1 - e_out) x [(1 + q_out)(1 + e_out) / sqrt(1 - e_out)]^(2/5) | Mardling & Aarseth 2001, as quoted in https://arxiv.org/abs/2306.09400 (original not opened) |

## Planets

| What | Value | Source |
|---|---|---|
| Planets per star, 0.5 to 10 AU | 1.6 (+0.72 / -0.89) | Cassan et al. 2012, https://arxiv.org/abs/1202.0903 |
| Small planets of M dwarfs (1-4 Earth radii, period under 200 days) | 2.5 +- 0.2 per star | Dressing & Charbonneau 2015, https://arxiv.org/abs/1501.01623 |
| Giant planets of sun-like stars | 10.5% (0.3-10 Jupiter masses, period 2-2000 days) | Cumming et al. 2008, https://arxiv.org/abs/0803.3357 |
| Giants against the star's mass and metals | share goes as mass x 10^(1.2 [Fe/H]); 3.5% at 0.5 Msun, 14% at 2 Msun | Johnson et al. 2010, https://arxiv.org/abs/1005.3084 |
| Hot Jupiters | 1.2 +- 0.38 % of sun-like stars | Wright et al. 2012, https://arxiv.org/abs/1205.2273 |
| Earth-size planets in the habitable zone | sun-like: 0.37 to 0.60 per star; M dwarfs: 0.16 to 0.24 | Bryson et al. 2021, https://arxiv.org/abs/2010.14812; Dressing & Charbonneau 2015 |
| The habitable zone: where a rocky planet keeps liquid water | 0.99 to 1.70 AU for the Sun (water loss, maximum greenhouse); for another star, times the square root of its power. Their fit to the star's temperature is not used | Kopparapu et al. 2013, https://arxiv.org/abs/1301.6674 (the older limits, Kasting et al. 1993: 0.95 to 1.67 AU) |
| Spacing of neighbours | most often 20 mutual Hill radii, 93% of pairs over 10, no period ratio under 1.2 | Weiss et al. 2018, https://arxiv.org/abs/1706.06204; Fabrycky et al. 2014, https://arxiv.org/abs/1202.6328 |
| Planets of double stars | around one star: a / a_binary < 0.464 - 0.380 mu - 0.631 e + 0.586 mu e + 0.150 e^2 - 0.198 mu e^2; around both: a / a_binary > 1.60 + 5.10 e - 2.22 e^2 + 4.12 mu - 4.27 e mu - 5.09 mu^2 + 4.61 e^2 mu^2 | Holman & Wiegert 1999, https://arxiv.org/abs/astro-ph/9809315 |
| Known planets | 6470 in the table `pscomppars` on 2026-11-15 | NASA Exoplanet Archive, https://exoplanetarchive.ipac.caltech.edu/TAP/sync?query=select+*+from+pscomppars&format=csv |

## Rare objects and other populations

| What | Value | Source |
|---|---|---|
| Neutron stars / stellar black holes | 1e9 / 1e8 (black holes: 1e7 to 1e9) | https://en.wikipedia.org/wiki/Milky_Way |
| Radio pulsars beamed at us | 30,000 +- 1,100 (the total, corrected for beaming: not verified) | Lorimer et al. 2006, https://arxiv.org/abs/astro-ph/0607640 |
| Magnetars | 31 known (how many are active in all: not verified) | https://en.wikipedia.org/wiki/Magnetar |
| X-ray binaries | about 350 catalogued (the whole population: not verified) | https://en.wikipedia.org/wiki/X-ray_binary; https://arxiv.org/abs/2302.02656 |
| Hypervelocity stars | about 1000 | Brown et al. 2007, via https://en.wikipedia.org/wiki/Hypervelocity_star |
| O stars | about 20,000 | https://en.wikipedia.org/wiki/O-type_star |
| Brown dwarfs | one for six stars near the Sun | https://en.wikipedia.org/wiki/Brown_dwarf |
| Free-floating planets | 21 (+23 / -13) per star, 0.33 to 6660 Earth masses | Sumi et al. 2023, https://arxiv.org/abs/2303.08280 |

## Stars that have died

| What | Value | Source |
|---|---|---|
| Life on the main sequence | 1e10 years x (mass in solar masses)^-2.5 | https://en.wikipedia.org/wiki/Main_sequence |
| What is left | under about 8 Msun a white dwarf, 9 to 25 a neutron star, over 25 a black hole (at the Sun's metals; not a clean cut, Sukhbold et al. 2016, https://arxiv.org/abs/1510.04643) | https://en.wikipedia.org/wiki/Supernova |
| Age of the thin disc | 8.8 +- 1.7 billion years | del Peloso et al. 2005, via https://en.wikipedia.org/wiki/Thin_disk |

## The ground of the Earth (`bodies.planetObjects.Earth`, `grids/earth_height`)

- **Heights**: ETOPO 2022 Global Relief Model, 60 arc-second, ice surface version
  (`ETOPO_2022_v1_60s_N90W180_surface.nc`, 478 MB), NOAA National Centers for Environmental
  Information, doi:10.25921/fd45-gt74. Land heights and sea depths in metres. NOAA says it is
  free for any scientific, personal or commercial use, excluding navigation. Made into the
  viewer's grid by `tools/make_terrain.py` at half the resolution (10800 x 5400, about 3.7 km a
  cell at the equator): the mean of two by two cells of the source.
- **Climate** (what grows where, `planetStyles.styles.earth` and the constants of `moisture` in
  `viewer/src/planet.wgsl`): not measured. A model of temperature (latitude, height) and rain
  (sea upwind, mountains in the way, belts of rising and sinking air), whose numbers were fitted
  by `tools/fit_climate.py` so that the colors of the land come out as in a cloudless picture of
  the Earth (the day cube map of the Gaia Sky data). The picture is only compared with; none of
  it is in the map.
- **Air**: the shares of red, green and blue scattered are those of Rayleigh scattering
  (wavelength to the minus fourth power, for 680, 550 and 440 nm); the thickness (0.27 for blue,
  straight down) is not verified against a source.

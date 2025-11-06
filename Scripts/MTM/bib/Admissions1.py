import matplotlib.pyplot as plt

SCHOOL_SYNONYMS = {
    # Texas A&M
    "texas a&m": "Texas A&M University",
    "a&m": "Texas A&M University",

    # University of Michigan
    "university of michigan": "University of Michigan",

    # Northwestern
    "northwestern university": "Northwestern University",

    # University of Chicago
    "university of chicago": "University of Chicago",
    "the university of chicago": "University of Chicago",

    # Vanderbilt
    "vanderbilt university": "Vanderbilt University",

    # Caltech
    "california institute of technology": "California Institute of Technology",
    "california institute of technology ": "California Institute of Technology",  # trailing space variant
    "caltech": "California Institute of Technology",

    # Rice
    "rice university": "Rice University",
    "rice chemistry": "Rice University",

    # Georgia Tech
    "georgia institute of technology": "Georgia Institute of Technology",
    "georgia tech": "Georgia Institute of Technology",

    # UC Santa Barbara
    "university of california santa barbara": "University of California, Santa Barbara",

    # UC San Diego
    "ucsd": "University of California San Diego",

    # EPFL
    "epfl": "EPFL",

    # University of Utah
    "university of utah": "University of Utah",

    # UW–Madison (normalize all forms)
    "uw-madison": "University of Wisconsin–Madison",
    "uw madison": "University of Wisconsin–Madison",
    "university of wisconsin - madison": "University of Wisconsin–Madison",
    "university of wisconsin madison": "University of Wisconsin–Madison",

    # Duke
    "duke university": "Duke University",

    # Stanford
    "stanford university": "Stanford University",

    # UT San Antonio
    "university of texas at san antonio": "University of Texas at San Antonio",

    # CU Boulder
    "cu boulder": "University of Colorado Boulder",
    "university of colorado boulder": "University of Colorado Boulder",

    # NYU
    "new york university": "New York University",
    "nyu": "New York University",

    # UCLA
    "ucla": "University of California, Los Angeles",

    # University of Washington (Seattle)
    "university of washington - seattle": "University of Washington (Seattle)",

    # Cornell
    "cornell university": "Cornell University",
    "cornell": "Cornell University",

    # Florida State
    "florida state university": "Florida State University",

    # Yale
    "yale": "Yale University",

    # UNC Chapel Hill
    "university of north carolina at chapel hill": "University of North Carolina at Chapel Hill",
    "unc-chapel hill": "University of North Carolina at Chapel Hill",
    "unc chapel hill": "University of North Carolina at Chapel Hill",

    # UC Berkeley
    "uc berkeley": "University of California, Berkeley",
    "university of california, berkeley": "University of California, Berkeley",

    # USC
    "university of southern california": "University of Southern California",
    "the university of southern california": "University of Southern California",

    # University of Minnesota Twin Cities
    "university of minnesota—twin cities": "University of Minnesota Twin Cities",
    "university of minnesota twin cities": "University of Minnesota Twin Cities",

    # UPenn
    "upenn": "University of Pennsylvania",

    # MIT
    "mit": "Massachusetts Institute of Technology",

    # Emory
    "emory university": "Emory University",

    # Nebraska–Lincoln
    "university of nebraska-lincoln": "University of Nebraska–Lincoln",

    # UIUC
    "university of illinois at urbana-champaign": "University of Illinois at Urbana-Champaign",

    # Penn State
    "penn state university": "Pennsylvania State University",
}

def normalize_school(name: str) -> str:
    key = (name or "").strip().lower()
    return SCHOOL_SYNONYMS.get(key, name.strip())

raw_list = [
    "Texas A&M", "University of Michigan", "Northwestern University", "University of Chicago",
    "Vanderbilt University", "The University of Chicago", "California Institute of Technology",
    "Rice University", "Georgia Institute of Technology", "Georgia Tech",
    "University of California Santa Barbara", "UCSD", "EPFL", "University of Utah",
    "UW-Madison", "Duke University", "Stanford University", "University of Texas at San Antonio",
    "Cu boulder", "New York University", "Texas A&M", "UCLA", "uw-madison",
    "Northwestern University", "University of Washington - Seattle", "Cornell University",
    "Caltech", "Cornell", "UW Madison", "NYU", "University of Wisconsin - Madison",
    "Florida State University", "Yale", "University of North Carolina at Chapel Hill",
    "UC Berkeley", "University of Southern California", "University of Minnesota—Twin Cities",
    "UPenn", "Stanford University", "Northwestern University", "MIT", "University of Utah",
    "UCLA", "University of North Carolina at Chapel Hill", "Cornell", "Emory University",
    "Cornell University", "California Institute of Technology ", "Georgia Tech",
    "UNC-Chapel Hill", "University of California, Berkeley", "Northwestern University",
    "University of Wisconsin Madison", "University of Wisconsin - Madison", "UCLA",
    "University of Nebraska-Lincoln", "UNC Chapel Hill", "Northwestern University",
    "University of Colorado Boulder", "Rice Chemistry", "University of Illinois at Urbana-Champaign ",
    "A&M ", "University of California, Berkeley", "UW Madison", "Duke University",
    "The University of Southern California", "Duke University", "Cornell University",
    "Penn State University", "University of Michigan", "University of California, Berkeley",
]

# Normalize and produce a canonical histogram
from collections import Counter
normalized = [normalize_school(x) for x in raw_list]
hist = Counter(normalized)
print(dict(hist))

institutions = dict(hist)

#
# institutions = {
#     "Northwestern Univ.": 13,
#     "UC Berkeley": 10,
#     "Cornell Univ.": 8,
#     "UCLA": 7,
#     "Texas A&M": 6,
#     "UW-Madison": 6,
#     "Stanford Univ.": 5,
#     "Georgia Tech": 4,
#     "Univ. of Michigan": 4,
#     "Caltech": 4,
#     "Univ. of North Carolina": 4,
#     "Univ. of Chicago": 4,
#     "Duke Univ.": 3,
#     "Rice Univ.": 3,
#     "Univ. of Washington": 2,
#     "Vanderbilt Univ.": 2,
#     "Univ. of Southern California": 2,
#     "Univ. of Utah": 2,
#     "UT San Antonio": 2,
#     "UPenn": 2,
#     "NYU": 2,
#     "Emory Univ.": 1,
#     "Florida State Univ.": 1,
#     "Univ. of Nebraska-Lincoln": 1,
#     "Univ. of Colorado Boulder": 1,
#     "UIUC": 1,
#     "UC Santa Barbara": 1,
#     "UCSD": 1,
#     "EPFL": 1,
#     "Penn State Univ.": 1
# }

school_colors = {
    "Northwestern University": "#4E2A84",
    "University of California, Berkeley": "#003262",
    "Cornell University": "#B31B1B",
    "University of California, Los Angeles": "#2774AE",
    "Texas A&M University": "#500000",
    "University of Wisconsin–Madison": "#C5050C",
    "Stanford University": "#8C1515",
    "Georgia Institute of Technology": "#B3A369",
    "University of Michigan": "#00274C",
    "California Institute of Technology": "#FF6C00",
    "University of North Carolina at Chapel Hill": "#7BAFD4",
    "University of Chicago": "#800000",
    "Duke University": "#001A57",
    "Rice University": "#00205B",
    "University of Washington (Seattle)": "#4B2E83",
    "Vanderbilt University": "#866D4B",
    "University of Southern California": "#990000",
    "University of Utah": "#CC0000",
    "University of Texas at San Antonio": "#F47321",
    "New York University": "#57068C",
    "University of Pennsylvania": "#011F5B",
    "Emory University": "#011E41",
    "Florida State University": "#782F40",
    "University of Nebraska–Lincoln": "#D00000",
    "University of Colorado Boulder": "#CFB87C",
    "University of Illinois at Urbana-Champaign": "#13294B",
    "University of California, Santa Barbara": "#003660",
    "University of California San Diego": "#00629B",
    "EPFL": "#E1302A",
    "Pennsylvania State University": "#001E44",
    "Yale University": "#0F4D92",
    "University of Minnesota Twin Cities": "#7A0019",
    "Massachusetts Institute of Technology": "#A31F34",
}

rankings = {}



# Sort for plotting
institutions_sorted = dict(sorted(institutions.items(), key=lambda item: item[1], reverse=True))
colors = [school_colors[school] for school in institutions_sorted.keys()]

# Print total number in different institutions
total_students = sum(institutions.values())
print(f"Total students choosing other institutions over UT Austin: {total_students}")

plt.figure(figsize=(10, 8))
plt.barh(list(institutions_sorted.keys()), list(institutions_sorted.values()), color=colors)
plt.xlabel("Number of Students")
plt.title("Institutions Chosen Over UT Austin")
plt.gca().invert_yaxis()
plt.tight_layout()
plt.show()
from scholarly import scholarly as sc

sq = sc.search_author("Jeanne Pemberton, Arizona")
print(sq)
print(next(sq))

exit()
my_id = 'hTpHqCYAAAAJ'

author = sc.search_author_id(my_id)
pubs=author.fill().publications
years = [p.bib['year'] for p in pubs]

#publist = [p.fill() for p in pubs]

#print(publist[0])
#print(len(publist))
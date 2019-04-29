import pickle
mediumCollege=pickle.load(open("mediumCollege.pickle","rb"))
mediumEmployer=pickle.load(open("mediumEmployer.pickle","rb"))
mediumLocation=pickle.load(open("mediumLocation.pickle","rb"))
print(mediumLocation)
print(mediumEmployer)
print(mediumCollege)
print(len(mediumCollege))


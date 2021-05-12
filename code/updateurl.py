#!/usr/env/bin python
"""A script to update the URLs of the NCMA data

"""
with open('ncma.txt') as f:
    data = json.load(f)
    newdata = []
    for rec in data["aaData"]:
        newrec = rec
        newrec[9] = '<a href=\"https://ncma.bigelow.org' + "/"  + newrec[0] +'">NCMA</a>'
        newdata.append(rec)
    data["aaData"] = newdata
    with open('ncma.json', 'w') as fout:
        json.dump(data, fout)

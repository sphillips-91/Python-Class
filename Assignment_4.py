import pandas as pd

data = {'Age': list(range(0, 101)), 'Fatality Rate': [20/1000000]*18 + [500/1000000]*32 + [6000/1000000]*15 + [90000/1000000]*36}
dt1 = pd.DataFrame(data)
dt2 = pd.read_csv("WorldDemographics.csv", index_col=0)
merged_dt = pd.merge(dt2, dt1, on="Age", how = "left")
merged_dt["Expected Deaths"] = merged_dt["#Alive"] * merged_dt["Fatality Rate"]
dt3 = merged_dt.groupby("PopulationID").agg({"#Alive": "sum", "Expected Deaths": "sum"})
dt3["Percent Died"] = (dt3["Expected Deaths"] / dt3["#Alive"]) * 100
dt3.to_csv("Assignment4_csv")
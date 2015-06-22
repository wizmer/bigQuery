frame = pd.DataFrame.from_csv("datasets/R_resolution.csv")
frame = frame.divide(frame.sum(axis=1),axis=0)
rgdtMeasured  = np.array(frame.columns.astype(float))
rgdtTheoretic = np.array(frame.index.astype(float))
rgdtF = frame.values[:-1,:-1]

frame = pd.DataFrame.from_csv("datasets/B_resolution.csv")
frame = frame.divide(frame.sum(axis=1),axis=0)
betaMeasured  = np.array(frame.columns.astype(float))
betaTheoretic = np.array(frame.index.astype(float))
betaF = frame.values[:-1,:-1]

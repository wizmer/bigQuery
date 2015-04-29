import bigQueryPlotting as bq

masks=[]
masks.append("notFirstTwo")
masks.append("notInSaaCut")
masks.append("zenithCut")
masks.append("runtypeCut")
masks.append("oneTRDTrack")
masks.append("goldenTRACKER")
masks.append("oneTrack")
masks.append("chargeOne")
masks.append("downGoing")
masks.append("betaNotCrazy")

theMask=bq.makeSelectionMask(masks)
print "{0:b}".format(theMask)

mass = "(Rfull * sqrt(1 - pow(BetaTOF,2)) / BetaTOF)"

whereClause="(Rfull > 0 && (selStatus&" + str(theMask)+ ")==" + str(theMask)+ " && " + mass + " > 0.8 && " + mass + " < 1.3 )"

isPhysicsTrigger=" (PhysBPatt >> 1)&"+str(int('11111',2))+ " as isPhysicsTrigger "
isTof="(JMembPatt>>5)&1 as isTof"
isEcal="(JMembPatt>>11)&1 as isEcal"

variables='{},{},{}'.format(isPhysicsTrigger,isTof,isEcal)

theCommand="bq query 'SELECT "+ variables + " from [Preselected.Preselected] WHERE " + whereClause + " GROUP BY isPhysicsTrigger, isTof, isEcal'"
print theCommand
print(bq.executeQuery(theCommand))
        
    
                                        

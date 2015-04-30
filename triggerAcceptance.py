import bigQueryPlotting as b

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


theTable="full_test.full_test"
b.setTable(theTable)

theMask=b.makeSelectionMask(masks)
print "{0:b}".format(theMask)

mass = "(Rfull * sqrt(1 - pow(BetaTOF,2)) / BetaTOF)"

nBins=20
firstBin=0
lastBin=20

#whereClause="(Rfull > 0 && (selStatus&" + str(theMask)+ ")==" + str(theMask)+ " && " + mass + " > 0.8 && " + mass + " < 1.3 )"
whereClause="(Rfull > 0 && (selStatus&" + str(theMask)+ ")==" + str(theMask)+ " )"

isPhysicsTrigger=" (PhysBPatt >> 1)&"+str(int('11111',2))+ " != 0 as isPhysicsTrigger "
isTof="(JMembPatt>>4)&1 as isTof"
isEcal="(JMembPatt>>11)&1 as isEcal"

variables='{} as binX, {},{},{},COUNT(1)'.format(b.bin(nBins,firstBin,lastBin,'Rfull'),isPhysicsTrigger,isTof,isEcal)

theCommand="""bq --format json query -n 10000 '
SELECT binX, IF(nTof + nEcal > 0,nPhysics*100/(nPhysics + nEcal*1000 + nTof*100),100), nPhysics, nEcal, nTof FROM (

    SELECT binX,SUM(IF(hashCode==221,CAST(f0_ AS INTEGER),0)) AS nPhysics,SUM(IF(hashCode==210,CAST(f0_ AS INTEGER),0)) AS nEcal, SUM(IF(hashCode==100,CAST(f0_ AS INTEGER),0)) AS nTof FROM

        (SELECT binX,isPhysicsTrigger,isEcal,isTof,(isPhysicsTrigger+10*IFNULL(CAST(isEcal AS INTEGER),INTEGER(2))+100*IFNULL(CAST(isTof AS INTEGER),INTEGER(2))) as hashCode,f0_  FROM

            (SELECT """ + variables  + """ FROM [""" + theTable + """] WHERE """ + whereClause + """ GROUP BY ROLLUP (binX,isPhysicsTrigger, isEcal, isTof) HAVING (isPhysicsTrigger==True && isTof IS NULL && isEcal IS NULL) || (isPhysicsTrigger==False && isEcal==True && isTof IS NULL) || (isPhysicsTrigger==False && isEcal==False && isTof==True) ORDER BY binX)

        ) GROUP BY binX )'"""

#print(b.executeQuery(theCommand))

b.histCustomCommand( nBins, firstBin, lastBin, theCommand)
    
#b.hist(100,0,1,'BetaTOF')
    
                                        

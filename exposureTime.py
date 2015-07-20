import bigQueryPlotting as b
import matplotlib.pyplot as plt

theTable="AMS.cutoffs"

queryOption=" --format json "
#queryOption=str()
globalOptions=str()

theCommand="""bq """ + globalOptions + """ query """ + queryOption + """'
    SELECT
      SUM(Lifetime) OVER (ORDER BY binX),
      cut as binX
    FROM (
      SELECT
        NTH_VALUE(Lifetime,1) OVER(ORDER BY cut) AS total,
        SUM(Lifetime) AS Lifetime,
        FLOOR(IGRF40pos) + 1 AS cut
      FROM """ + theTable + """
      WHERE (
        Lifetime > 0.5 &&
        ABS(ThetaM) > 0.5
      )
      GROUP BY
        ROLLUP (cut)
      ORDER BY
        cut )
    WHERE
      cut > 0
'"""




#data=b.histCustomCommand( nBins, firstBin, lastBin, theCommand)
df=b.histCustomCommand(theCommand)
h=b.Hist( df, 1, 0, 50 )
h.plot()
plt.show()

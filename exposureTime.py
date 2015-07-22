import bigQueryPlotting as b
import matplotlib.pyplot as plt

theTable="AMS.cutoffs"
queryOption=str()
globalOptions=str()

queryOption=" --format json "

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
      JOIN AMS.timeInSecInData
      ON (AMS.timeInSecInData.JMDCTimeInSec = """ + theTable + """.JMDCTime)
      WHERE (
         goodSecond == 1
      )
      GROUP BY
        ROLLUP (cut)
      ORDER BY
        cut )
    WHERE
      cut > 0
'"""

df=b.histCustomCommand(theCommand)
h=b.Hist( df, 1, 0, 50 )
h.plot()
plt.show()

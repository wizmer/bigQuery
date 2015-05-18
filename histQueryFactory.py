class HistQueryFactory(object):
    """
    h = HistQueryFactory()
    h.add_variable("R",0,10,40)
    h.add_variable("BetaTOF",0,1,10)
    h.add_condition("Rcutoff < R")
    h.add_condition("Latitude > 0.8")
    h.add_condition("Latitude < 0.9")
    print str(h)
    """    
    def __init__(self, table="full_test.full_test", bins='mid'):
        self.table = table
        self.bins = bins
        self.variables = {}
        self.varorder = []
        self.extra_conditions = []
    
    def add_variable(self, var, min=0, max=1, nbins=100):
        self.variables[var] = [min, max, nbins]
        self.varorder.append(var)
    
    def add_condition(self, condition):
        self.extra_conditions.append(condition)
    
    def __str__(self):
        select, where = [], []
        binexpr = "STRING({mi:f} + ({rng:f})*FLOOR({N:f}*({var}-({mi:f}))/{rng:f})/{N:f}) as {var}_BIN"
        for v in self.varorder:
            mi,ma,nb = self.variables[v]
            select.append(binexpr.format(var=v, mi=mi, N=nb, rng=ma-mi))
            where.append("{0} >= {1} AND {0} < {2}".format(v,mi,ma))
        
        select.append("COUNT(1) as count")
        where += self.extra_conditions
        
        query = \
        "SELECT\n   {select}\n"\
        "FROM\n   {from}\n"\
        "WHERE\n   {where}\n"\
        "GROUP BY {bins}\n"\
        "ORDER BY {bins}\n"
        
        dic = {
            "select": ",\n   ".join(select),
            "from": self.table,
            "where":" AND\n   ".join(where),
            "bins":",".join([v+"_BIN" for v in self.varorder])
        }
        return  query.format(**dic)    

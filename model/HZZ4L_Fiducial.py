from HiggsAnalysis.CombinedLimit.PhysicsModel import *
# from HiggsAnalysis.CombinedLimit.SMHiggsBuilder import SMHiggsBuilder

class InclusiveFiducial( PhysicsModel ):
    ''' Model used to unfold differential distributions '''

    def __init__(self):
        PhysicsModel.__init__(self)
        self.Range=[0.,10]
        self.fracRange=[0.,0.5]                
        self.mHRange=[20,1000]
        self.mass=125.0
        self.debug=1

    def setPhysicsOptions(self,physOptions):
        if self.debug>0:print "Setting PhysicsModel Options"
        for po in physOptions:
            if po.startswith("range="):
                self.Range=po.replace("range=","").split(",")
                if len(self.Range)!=2:
                    raise RunTimeError, "Range require minimal and maximal values: range=min,max"
                if self.debug>0:print "new Range is ", self.Range
            if po.startswith("nBin="):
                self.nBin=int(po.replace("nBin=",""))
                if self.debug>0:print "new n. of bins is ",self.nBin
            if po.startswith("higgsMassRange="):
                if self.debug>0: print "setting higgs mass range floating:",po.replace("higgsMassRange=","").split(",")
                self.mHRange=po.replace("higgsMassRange=","").split(",")
                #checks
                if len(self.mHRange) != 2:
                    raise RuntimeError, "Higgs mass range definition requires two extrema"
                elif float(self.mHRange[0]) >= float(self.mHRange[1]):
                    raise RuntimeError, "Extrema for Higgs mass range defined with inverterd order. Second must be larger the first"
            if po.startswith("mass="):
                self.mass=float( po.replace('mass=','') )

            #verbose
            if po.startswith("verbose"):
                self.debug = 1

    def doParametersOfInterest(self):
        POIs=""
        if self.debug>0:print "Setting pois"
                
        self.modelBuilder.doVar("r[1,%s,%s]" % (self.Range[0], self.Range[1]))
        self.modelBuilder.doVar("frac4e[0.25,%s,%s]" % (self.fracRange[0], self.fracRange[1]))
        self.modelBuilder.doVar("frac4mu[0.25,%s,%s]" % (self.fracRange[0], self.fracRange[1]))
                
        POIs+="r,"
        POIs+="frac4e,"
        POIs+="frac4mu"
                                        
        poiNames=[]
        if self.modelBuilder.out.var("MH"):
            if len(self.mHRange) == 2:
                print 'MH will be left floating within', self.mHRange[0], 'and', self.mHRange[1]
                self.modelBuilder.out.var("MH").setRange(float(self.mHRange[0]),float(self.mHRange[1]))
                self.modelBuilder.out.var("MH").setConstant(False)
                poiNames += [ 'MH' ]
            else:
                print 'MH will be assumed to be', self.mass
                self.modelBuilder.out.var("MH").removeRange()
                self.modelBuilder.out.var("MH").setVal(self.mass)
        else:
            if len(self.mHRange) == 2:
                print 'MH will be left floating within', self.mHRange[0], 'and', self.mHRange[1]
                self.modelBuilder.doVar("MH[%s,%s]" % (self.mHRange[0],self.mHRange[1]))
                poiNames += [ 'MH' ]
            else:
                print 'MH (not there before) will be assumed to be', self.mass
                self.modelBuilder.doVar("MH[%g]" % self.mass)
        for poi in poiNames:
                POIs += ",%s"%poi
        self.modelBuilder.doSet("POI",POIs)
        print "set up pois"
        self.setup()

    def setup(self):
        self.modelBuilder.factory_('expr::Sigma_trueH4e("@0*@1", r, frac4e)')
        self.modelBuilder.factory_('expr::Sigma_trueH4mu("@0*@1", r, frac4mu)')
        self.modelBuilder.factory_('expr::Sigma_trueH2e2mu("@0*(1-@1-@2)", r, frac4e, frac4mu)')
        self.modelBuilder.factory_('expr::Sigma_trueZ4e("@0*@1", r, frac4e)')
        self.modelBuilder.factory_('expr::Sigma_trueZ4mu("@0*@1", r, frac4mu)')
        self.modelBuilder.factory_('expr::Sigma_trueZ2e2mu("@0*(1-@1-@2)", r, frac4e, frac4mu)')


    def getYieldScale(self,bin,process):
        if not self.DC.isSignal[process]: return 1
        if process in [ "trueH4e","trueH4mu","trueH2e2mu","trueZ4e","trueZ4mu","trueZ2e2mu"]: 
            return 'Sigma_'+process
        else : return 1


class DifferentialFiducial( PhysicsModel ):
    ''' Model used to unfold differential distributions '''

    def __init__(self):
        PhysicsModel.__init__(self)
        self.Range=[0.,10]
        self.fracRange=[0.,0.5]
        self.nBin=4
        self.mHRange=[20.0,200.0]
        self.debug=1
        self.mass=125.0

    def setPhysicsOptions(self,physOptions):
        if self.debug>0:print "Setting PhysicsModel Options"
        for po in physOptions:
            if po.startswith("range="):
                self.Range=po.replace("range=","").split(",")
                if len(self.Range)!=2:
                    raise RunTimeError, "Range require minimal and maximal values: range=min,max"
                if self.debug>0:print "new Range is ", self.Range
            if po.startswith("nBin="):
                self.nBin=int(po.replace("nBin=",""))
                if self.debug>0:print "new n. of bins is ",self.nBin
            if po.startswith("higgsMassRange="):
                if self.debug>0: print "setting higgs mass range floating:",po.replace("higgsMassRange=","").split(",")
                self.mHRange=po.replace("higgsMassRange=","").split(",")
                #checks
                if len(self.mHRange) != 2:
                    raise RuntimeError, "Higgs mass range definition requires two extrema"
                elif float(self.mHRange[0]) >= float(self.mHRange[1]):
                    raise RuntimeError, "Extrema for Higgs mass range defined with inverterd order. Second must be larger the first"
            if po.startswith("mass="):
                self.mass=float( po.replace('mass=','') )
            #verbose
            if po.startswith("verbose"):
                self.debug = 1

    def doParametersOfInterest(self):
        POIs=""
        if self.debug>0:print "Setting pois"
        for iBin in range(0,self.nBin):
            fracSM4e = self.modelBuilder.out.var("fracSM4eBin%d" % (iBin)).getVal()
            fracSM4mu = self.modelBuilder.out.var("fracSM4muBin%d" % (iBin)).getVal()

            if self.modelBuilder.out.var("rBin%d" % (iBin)):
                self.modelBuilder.out.var("rBin%d" % (iBin)).setRange(self.Range[0], self.Range[1])
                self.modelBuilder.out.var("rBin%d" % (iBin)).setConstant(False)
            else :
                self.modelBuilder.doVar("rBin%d[1, %s,%s]" % (iBin, self.Range[0],self.Range[1]))
            ''''
            if self.modelBuilder.out.var("fracSM4eBin%d" % (iBin)):
                self.modelBuilder.out.var("fracSM4eBin%d" % (iBin)).setRange(self.fracRange[0], self.fracRange[1])
                self.modelBuilder.out.var("fracSM4eBin%d" % (iBin)).setConstant(False)
            else :
                self.modelBuilder.doVar("fracSM4eBin%d[0.25, %s,%s]" % (iBin, self.fracRange[0],self.fracRange[1]))

            if self.modelBuilder.out.var("fracSM4muBin%d" % (iBin)):
                self.modelBuilder.out.var("fracSM4muBin%d" % (iBin)).setRange(self.fracRange[0], self.fracRange[1])
                self.modelBuilder.out.var("fracSM4muBin%d" % (iBin)).setConstant(False)
            else :
                self.modelBuilder.doVar("fracSM4muBin%d[0.25, %s,%s]" % (iBin, self.fracRange[0],self.fracRange[1]))
            '''
            if iBin>=0:
                POIs+="rBin%d,"%iBin
                #POIs+="fracSM4eBin%d,"%iBin
                #POIs+="fracSM4muBin%d,"%iBin
                if self.debug>0:print "Added Bin%d to the POIs"%iBin
        poiNames=[]
        if self.modelBuilder.out.var("MH"):
            if len(self.mHRange) == 2:
                print 'MH will be left floating within', self.mHRange[0], 'and', self.mHRange[1]
                self.modelBuilder.out.var("MH").setRange(float(self.mHRange[0]),float(self.mHRange[1]))
                self.modelBuilder.out.var("MH").setConstant(False)
                poiNames += [ 'MH' ]
            else:
                print 'MH will be assumed to be', self.mass
                self.modelBuilder.out.var("MH").removeRange()
                self.modelBuilder.out.var("MH").setVal(self.mass)
        else:
            if len(self.mHRange) == 2:
                print 'MH will be left floating within', self.mHRange[0], 'and', self.mHRange[1]
                self.modelBuilder.doVar("MH[%s,%s]" % (self.mHRange[0],self.mHRange[1]))
                poiNames += [ 'MH' ]
            else:
                print 'MH (not there before) will be assumed to be', self.mass
                self.modelBuilder.doVar("MH[%g]" % self.mass)
        for poi in poiNames:
            POIs += "%s,"%poi
        POIs = POIs[:-1] # remove last comma
        self.modelBuilder.doSet("POI",POIs)
        self.setup()

    def setup(self):        
        for iBin in range(0,self.nBin):
            #these define the signal strenghts per production mode and should be implemented in createXSworkspace
            self.modelBuilder.factory_('expr::r_trueH2e2muBin%d("@0*(1-@1-@2)", rBin%d, fracSM4eBin%d, fracSM4muBin%d)' % (iBin,iBin,iBin,iBin))
            self.modelBuilder.factory_('expr::r_trueH4eBin%d("@0*@1", rBin%d, fracSM4eBin%d)'% (iBin,iBin,iBin))
            self.modelBuilder.factory_('expr::r_trueH4muBin%d("@0*@1", rBin%d, fracSM4muBin%d)'% (iBin,iBin,iBin))

    # def getYieldScale(self,bin,process):
    #     if not self.DC.isSignal[process]: return 1
    #     name = "fiducial_%s" % process
             
    #     self.modelBuilder.factory_('expr::%s("@0", r_%s)' % (name, process))
    #     if process in [ "trueH2e2muBin0","trueH4eBin0","trueH4muBin0","trueH2e2muBin1","trueH4eBin1","trueH4muBin1","trueH2e2muBin2","trueH4eBin2","trueH4muBin2","trueH2e2muBin3","trueH4eBin3","trueH4muBin3"]: 
    #         return name
    #     else :
    #         return 1
    
    def getYieldScale(self,bin,process):
        if not self.DC.isSignal[process]: return 1
        Processes = []
        ## Inclusive signal strength has to be 1.0
        ## i.e. sum of per fs mu has to be 1.0
        for iBin in range(0,self.nBin):
            for channel in ['4e', '4mu', '2e2mu']:
                Processes += ['trueH'+channel+'Bin'+str(iBin)]
        if process in Processes: return 'r_'+process
        else: return 1

class InclusiveFiducialV2( PhysicsModel ):
    ''' Model used to unfold differential distributions '''

    def __init__(self):
        PhysicsModel.__init__(self)
        self.SigmaRange=[0.,100]
        self.MHRange=[20,1000]
        self.defaultMH = 125.0
        self.debug=1

    def setPhysicsOptions(self,physOptions):
        if self.debug>0:print "Setting PhysicsModel Options"
        for po in physOptions:
            if po.startswith("range="):
                self.SigmaRange=po.replace("range=","").split(",")
                if len(self.SigmaRange)!=2:
                    raise RunTimeError, "SigmaRange require minimal and maximal values: range=min,max"
                if self.debug>0:print "New SigmaRange is ", self.SigmaRange
            if po.startswith("higgsMassRange="):
                if self.debug>0: print "setting MHRange floating:",po.replace("higgsMassRange=","").split(",")
                self.MHRange=po.replace("higgsMassRange=","").split(",")
                #checks
                if len(self.MHRange) != 2:
                    raise RuntimeError, "MHRange definition requires two extrema: higgsMassRange=min,max"
                elif float(self.MHRange[0]) >= float(self.MHRange[1]):
                    raise RuntimeError, "Extrema for MH defined with inverterd order. Second must be larger the first"
            if po.startswith("mass="):
                self.defaultMH=float( po.replace('mass=','') )
            #verbose
            if po.startswith("verbose"):
                self.debug = 1

    def doParametersOfInterest(self):
        POIs=""
        if self.debug>0:print "Setting POIs"


        if self.modelBuilder.out.var("Sigma4e"):
            self.modelBuilder.out.var("Sigma4e").setRange(self.SigmaRange[0], self.SigmaRange[1])
        else:
            self.modelBuilder.doVar("Sigma4e[1,%s,%s]" % (self.SigmaRange[0], self.SigmaRange[1]))
        if self.modelBuilder.out.var("Sigma4mu"):
            self.modelBuilder.out.var("Sigma4mu").setRange(self.SigmaRange[0], self.SigmaRange[1])
        else:
            self.modelBuilder.doVar("Sigma4mu[1,%s,%s]" % (self.SigmaRange[0], self.SigmaRange[1]))
        if self.modelBuilder.out.var("Sigma2e2mu"):
            self.modelBuilder.out.var("Sigma2e2mu").setRange(self.SigmaRange[0], self.SigmaRange[1])
        else:
            self.modelBuilder.doVar("Sigma2e2mu[1,%s,%s]" % (self.SigmaRange[0], self.SigmaRange[1]))

        POIs+="Sigma4e,"
        POIs+="Sigma4mu,"
        POIs+="Sigma2e2mu"

        poiNames=[]
        if self.modelBuilder.out.var("MH"):
            if len(self.MHRange) == 2:
                print 'MH will be left floating within', self.MHRange[0], 'and', self.MHRange[1]
                self.modelBuilder.out.var("MH").setRange(float(self.MHRange[0]),float(self.MHRange[1]))
                self.modelBuilder.out.var("MH").setVal(self.defaultMH)
                self.modelBuilder.out.var("MH").setConstant(False)
                poiNames += [ 'MH' ]
            else:
                print 'MH will be assumed to be', self.defaultMH
                self.modelBuilder.out.var("MH").removeRange()
                self.modelBuilder.out.var("MH").setVal(self.defaultMH)
                self.modelBuilder.out.var("MH").setConstant(True)
        else:
            if len(self.MHRange) == 2:
                print 'MH will be left floating within', self.MHRange[0], 'and', self.MHRange[1]
                self.modelBuilder.doVar("MH[%s,%s]" % (self.MHRange[0],self.MHRange[1]))
                self.modelBuilder.out.var("MH").setVal(self.defaultMH)
                poiNames += [ 'MH' ]
            else:
                print 'MH (not there before) will be assumed to be', self.defaultMH
                self.modelBuilder.doVar("MH[%g]" % self.defaultMH)

        for poi in poiNames:
            POIs += ",%s"%poi

        self.modelBuilder.doSet("POI",POIs)
        print "set up pois"
        self.setup()

    def setup(self):
        self.modelBuilder.factory_('expr::Sigma_trueH4e("@0", Sigma4e)')
        self.modelBuilder.factory_('expr::Sigma_trueH4mu("@0", Sigma4mu)')
        self.modelBuilder.factory_('expr::Sigma_trueH2e2mu("@0", Sigma2e2mu)')
        self.modelBuilder.factory_('expr::Sigma_trueZ4e("@0", Sigma4e)')
        self.modelBuilder.factory_('expr::Sigma_trueZ4mu("@0", Sigma4mu)')
        self.modelBuilder.factory_('expr::Sigma_trueZ2e2mu("@0", Sigma2e2mu)')

    def getYieldScale(self,bin,process):
        if not self.DC.isSignal[process]: return 1
        if process in [ "trueH4e","trueH4mu","trueH2e2mu","trueZ4e","trueZ4mu","trueZ2e2mu"]:
            return 'Sigma_'+process
        else:
            return 1


class DifferentialFiducialV2( PhysicsModel ):
    ''' Model used to unfold differential distributions '''

    def __init__(self):
        PhysicsModel.__init__(self)
        self.Range=[0.,10]
        self.fracRange=[0.,0.5]
        self.nBin=7
        self.mHRange=[]
        self.debug=1
        self.mass=0

    def setPhysicsOptions(self,physOptions):
        if self.debug>0:print "Setting PhysicsModel Options"
        for po in physOptions:
            if po.startswith("range="):
                self.Range=po.replace("range=","").split(":")
                if len(self.Range)!=2:
                    raise RunTimeError, "Range require minimal and maximal values: range=min:max"
                if self.debug>0:print "new Range is ", self.Range
            if po.startswith("nBin="):
                self.nBin=int(po.replace("nBin=",""))
                if self.debug>0:print "new n. of bins is ",self.nBin
            if po.startswith("higgsMassRange="):
                if self.debug>0: print "setting higgs mass range floating:",po.replace("higgsMassRange=","").split(":")
                self.mHRange=po.replace("higgsMassRange=","").split(",")
                #checks
                if len(self.mHRange) != 2:
                    raise RuntimeError, "Higgs mass range definition requires two extrema"
                elif float(self.mHRange[0]) >= float(self.mHRange[1]):
                    raise RuntimeError, "Extrema for Higgs mass range defined with inverterd order. Second must be larger the first"
            if po.startswith("mass="):
                self.mass=float( po.replace('mass=','') )

        #verbose
        if po.startswith("verbose"):
            self.debug = 1

    def doParametersOfInterest(self):
        POIs=""
        if self.debug>0:print "Setting pois"
        print "nBins:",self.nBin
        for iBin in range(0,self.nBin):
            print "bin",iBin    
            if self.modelBuilder.out.var("r2e2muBin%d" % (iBin)):  
                self.modelBuilder.out.var("r2e2muBin%d" % (iBin)).setRange(self.Range[0], self.Range[1])
                self.modelBuilder.out.var("r2e2muBin%d" % (iBin)).setConstant(False)
            else : 
                self.modelBuilder.doVar("r2e2muBin%d[1, %s,%s]" % (iBin, self.Range[0],self.Range[1]))

            if self.modelBuilder.out.var("r4muBin%d" % (iBin)):
                self.modelBuilder.out.var("r4muBin%d" % (iBin)).setRange(self.Range[0], self.Range[1])
                self.modelBuilder.out.var("r4muBin%d" % (iBin)).setConstant(False)
            else :
                self.modelBuilder.doVar("r4muBin%d[1, %s,%s]" % (iBin, self.Range[0],self.Range[1]))

            if self.modelBuilder.out.var("r4eBin%d" % (iBin)):
                self.modelBuilder.out.var("r4eBin%d" % (iBin)).setRange(self.Range[0], self.Range[1])
                self.modelBuilder.out.var("r4eBin%d" % (iBin)).setConstant(False)
            else :
                self.modelBuilder.doVar("r4eBin%d[1, %s,%s]" % (iBin, self.Range[0],self.Range[1]))

            if iBin>-1:
                POIs+="r2e2muBin%d,"%iBin
                POIs+="r4eBin%d,"%iBin
                POIs+="r4muBin%d,"%iBin
                if self.debug>0:print "Added Bin%d to the POIs"%iBin

        poiNames=[]
        if self.modelBuilder.out.var("MH"):
            if len(self.mHRange) == 2:
                print 'MH will be left floating within', self.mHRange[0], 'and', self.mHRange[1]
                self.modelBuilder.out.var("MH").setRange(float(self.mHRange[0]),float(self.mHRange[1]))
                self.modelBuilder.out.var("MH").setConstant(False)
                poiNames += [ 'MH' ]
            else:
                print 'MH will be assumed to be', self.mass
                self.modelBuilder.out.var("MH").removeRange()
                self.modelBuilder.out.var("MH").setVal(self.mass)
        else:
            if len(self.mHRange) == 2:
                print 'MH will be left floating within', self.mHRange[0], 'and', self.mHRange[1]
                self.modelBuilder.doVar("MH[%s,%s]" % (self.mHRange[0],self.mHRange[1]))
                poiNames += [ 'MH' ]
            else:
                print 'MH (not there before) will be assumed to be', self.mass
                self.modelBuilder.doVar("MH[%g]" % self.mass)
        for poi in poiNames:
            POIs += "%s,"%poi
        POIs = POIs[:-1] # remove last comma
        self.modelBuilder.doSet("POI",POIs)
        self.setup()

    def setup(self):
        for iBin in range(0,self.nBin):
            self.modelBuilder.factory_('expr::Sigma_trueH4eBin%d("@0", r4eBin%d)'% (iBin,iBin))
            self.modelBuilder.factory_('expr::Sigma_trueH4muBin%d("@0", r4muBin%d)'% (iBin,iBin))
            self.modelBuilder.factory_('expr::Sigma_trueH2e2muBin%d("@0", r2e2muBin%d)' % (iBin,iBin))
            self.modelBuilder.factory_('expr::Sigma_trueZ4eBin%d("@0", r4eBin%d)'% (iBin,iBin))
            self.modelBuilder.factory_('expr::Sigma_trueZ4muBin%d("@0", r4muBin%d)'% (iBin,iBin))
            self.modelBuilder.factory_('expr::Sigma_trueZ2e2muBin%d("@0", r2e2muBin%d)' % (iBin,iBin))

    def getYieldScale(self,bin,process):
        if not self.DC.isSignal[process]: return 1
        Processes = []
        for Boson in ['H', 'Z']:
            for iBin in range(0,self.nBin):
                for channel in ['4e', '4mu', '2e2mu']:
                    Processes += ['true'+Boson+channel+'Bin'+str(iBin)]
        if process in Processes: return 'Sigma_'+process
        else: return 1


class InclusiveFiducialV3( PhysicsModel ):
    ''' Fiducial cross-section model for both H->4l and Z->4l'''
# Define SM fractions, and set them as constants:
#  fracSM4e: const 
#  fracSM4mu: const
#  fracSM2e2mu: = 1 - fracSM4e - fracSM4mu
#
# Define parameters:
#  K1 : range [0, maxK1], where maxK1 =  1/fracSM4e
#  K2 : range [0, maxK2], where maxK2 = (1-fracSM4e)/fracSM4mu
#
# Define fid-xsecs of the total and 4e/4mu/2e2mu final states:
#  Sigma: total
#  Sigma4e: for 4e
#  Sigma4mu: for 4mu
#  Sigma2e2mu: for 2e2mu
#
# Express them in terms of Sigma total, SM fractions, and the K1, K2 parameters:
#  Sigma4e    = Sigma *    fracSM4e*K1
#  Sigma4mu   = Sigma * (1-fracSM4e*K1) *      K2 * fracSM4mu/(1-fracSM4e) 
#  Sigma2e2mu = Sigma * (1-fracSM4e*K1) * [1 - K2 * fracSM4mu/(1-fracSM4e)] 
#
# In this way, we have:
#                          Sigma4e           Sigma4mu                       Sigma2e2mu
# SM fractions:
#   K1=1,     K2=1     :   Sigma*fracSM4e    Sigma*fracSM4mu                Sigma*(1-fracSM4e-fracSM4mu)
# Extream cases:
#   K1=maxK1, K2=any   :   Sigma             0                              0
#   K1=0,     K2=maxK2 :   0                 Sigma                          0
#   K1=0,     K2=0     :   0                 0                              Sigma
# Other examples:
#   K1=1,     K2=maxK2 :   Sigma*fracSM4e    Sigma*(1-fracSM4e)             0
#   K1=1,     K2=0     :   Sigma*fracSM4e    0                              Sigma*(1-fracSM4e)
#   K1=0,     K2=1     :   0                 Sigma*fracSM4mu/(1-fracSM4e)   Sigma*[1-fracSM4mu/(1-fracSM4e)])

    def __init__(self):
        PhysicsModel.__init__(self)
        self.SigmaRange=[0.,100]
        self.MHRange=[20,1000]
        self.defaultMH = 125.0
        self.debug=1

    def setPhysicsOptions(self,physOptions):
        if self.debug>0:print "Setting PhysicsModel Options"
        for po in physOptions:
            if po.startswith("range="):
                self.SigmaRange=po.replace("range=","").split(",")
                if len(self.SigmaRange)!=2:
                    raise RunTimeError, "SigmaRange require minimal and maximal values: range=min,max"
                if self.debug>0:print "New SigmaRange is ", self.SigmaRange
            if po.startswith("higgsMassRange="):
                if self.debug>0: print "setting MHRange floating:",po.replace("higgsMassRange=","").split(",")
                self.MHRange=po.replace("higgsMassRange=","").split(",")
                #checks
                if len(self.MHRange) != 2:
                    raise RuntimeError, "MHRange definition requires two extrema: higgsMassRange=min,max"
                elif float(self.MHRange[0]) >= float(self.MHRange[1]):
                    raise RuntimeError, "Extrema for MH defined with inverterd order. Second must be larger the first"
            if po.startswith("mass="):
                self.defaultMH=float( po.replace('mass=','') )
            #verbose
            if po.startswith("verbose"):
                self.debug = 1

    def doParametersOfInterest(self):
        POIs=""
        if self.debug>0:print "Setting pois"
        
        # get values from the workspace
        fracSM4e = self.modelBuilder.out.var("fracSM4e").getVal()
        fracSM4mu = self.modelBuilder.out.var("fracSM4mu").getVal()

        # define parameters
        if self.modelBuilder.out.var("Sigma"):
            self.modelBuilder.out.var("Sigma").setRange(self.SigmaRange[0], self.SigmaRange[1])
        else:
            self.modelBuilder.doVar("Sigma[1,%s,%s]" % (self.SigmaRange[0], self.SigmaRange[1]))
        if self.modelBuilder.out.var("K1"):
            self.modelBuilder.out.var("K1").setRange(0.0, 1.0/fracSM4e)
        else:
            self.modelBuilder.doVar("K1[1,%s,%s]" % (0.0,  1.0/fracSM4e))
        if self.modelBuilder.out.var("K2"):
            self.modelBuilder.out.var("K2").setRange(0.0, (1.0-fracSM4e)/fracSM4mu)
        else:
            self.modelBuilder.doVar("K2[1,%s,%s]" % (0.0, (1.0-fracSM4e)/fracSM4mu))

        POIs+="Sigma,"
        ## Possible remove these, keep BR constant in the fit.
        ## Would make sense as they are always at boundary.
        POIs+="K1,"
        POIs+="K2"

        poiNames=[]
        if self.modelBuilder.out.var("MH"):
            if len(self.MHRange) == 2:
                print 'MH will be left floating within', self.MHRange[0], 'and', self.MHRange[1]
                self.modelBuilder.out.var("MH").setRange(float(self.MHRange[0]),float(self.MHRange[1]))
                self.modelBuilder.out.var("MH").setConstant(False)
                poiNames += [ 'MH' ]
            else:
                print 'MH will be assumed to be', self.defaultMH
                self.modelBuilder.out.var("MH").removeRange()
                self.modelBuilder.out.var("MH").setVal(self.defaultMH)
        else:
            if len(self.MHRange) == 2:
                print 'MH will be left floating within', self.MHRange[0], 'and', self.MHRange[1]
                self.modelBuilder.doVar("MH[%s,%s]" % (self.MHRange[0],self.MHRange[1]))
                poiNames += [ 'MH' ]
            else:
                print 'MH (not there before) will be assumed to be', self.defaultMH
                self.modelBuilder.doVar("MH[%g]" % self.defaultMH)
        for poi in poiNames:
            POIs += ",%s"%poi
        self.modelBuilder.doSet("POI",POIs)
        print "set up POIs"
        self.setup()

    def setup(self):
        self.modelBuilder.factory_('expr::Sigma_trueH4e("@0*@1*@2", Sigma, fracSM4e, K1)')
        self.modelBuilder.factory_('expr::Sigma_trueH4mu("@0*(1.0-@1*@2)*@3*@4/(1.0-@1)", Sigma, fracSM4e, K1, K2, fracSM4mu)')
        self.modelBuilder.factory_('expr::Sigma_trueH2e2mu("@0*(1.0-@1*@2)*(1.0-@3*@4/(1.0-@1))", Sigma, fracSM4e, K1, K2, fracSM4mu)')
        self.modelBuilder.factory_('expr::Sigma_trueZ4e("@0*@1*@2", Sigma, fracSM4e, K1)')
        self.modelBuilder.factory_('expr::Sigma_trueZ4mu("@0*(1.0-@1*@2)*@3*@4/(1.0-@1)", Sigma, fracSM4e, K1, K2, fracSM4mu)')
        self.modelBuilder.factory_('expr::Sigma_trueZ2e2mu("@0*(1.0-@1*@2)*(1.0-@3*@4/(1.0-@1))", Sigma, fracSM4e, K1, K2, fracSM4mu)')

    def getYieldScale(self,bin,process):
        if not self.DC.isSignal[process]: return 1
        if process in [ "trueH4e","trueH4mu","trueH2e2mu","trueZ4e","trueZ4mu","trueZ2e2mu"]:
            return 'Sigma_'+process
        else:
            return 1



class DifferentialFiducialV3( PhysicsModel ):
    ''' Model used to unfold differential distributions for Fiducial cross-section model for both H->4l and Z->4l'''

    def __init__(self):
        PhysicsModel.__init__(self)
        self.SigmaRange=[0.,100]
        self.nBin=4
        self.MHRange=[20.0,200.0]
        self.defautMH=125.0
        self.debug=1

    def setPhysicsOptions(self,physOptions):
        if self.debug>0:print "Setting PhysicsModel Options"
        for po in physOptions:
            if po.startswith("range="):
                self.SigmaRange=po.replace("range=","").split(",")
                if len(self.SigmaRange)!=2:
                    raise RunTimeError, "SigmaRange require minimal and maximal values: range=min,max"
                if self.debug>0:print "New SigmaRange is ", self.SigmaRange
            if po.startswith("higgsMassRange="):
                if self.debug>0: print "setting MHRange floating:",po.replace("higgsMassRange=","").split(",")
                self.MHRange=po.replace("higgsMassRange=","").split(",")
                #checks
                if len(self.MHRange) != 2:
                    raise RuntimeError, "MHRange definition requires two extrema: higgsMassRange=min,max"
                elif float(self.MHRange[0]) >= float(self.MHRange[1]):
                    raise RuntimeError, "Extrema for MH defined with inverterd order. Second must be larger the first"
            if po.startswith("mass="):
                self.defaultMH=float( po.replace('mass=','') )
            if po.startswith("nBin="):
                self.nBin=int(po.replace("nBin=",""))
                if self.debug>0:print "new n. of bins is ",self.nBin
            #verbose
            if po.startswith("verbose"):
                self.debug = 1


    def doParametersOfInterest(self):
        POIs=""
        if self.debug>0:print "Setting pois"

        for iBin in range(0,self.nBin):
            # get values from the workspace
            fracSM4e = self.modelBuilder.out.var("fracSM4eBin%d" % (iBin)).getVal()
            fracSM4mu = self.modelBuilder.out.var("fracSM4muBin%d" % (iBin)).getVal()

            if self.modelBuilder.out.var("SigmaBin%d" % (iBin)):
                self.modelBuilder.out.var("SigmaBin%d" % (iBin)).setRange(self.SigmaRange[0], self.SigmaRange[1])
                self.modelBuilder.out.var("SigmaBin%d" % (iBin)).setConstant(False)
            else :
                self.modelBuilder.doVar("SigmaBin%d[1, %s,%s]" % (iBi, self.SigmaRange[0],self.SigmaRange[1]))

            if self.modelBuilder.out.var("K1Bin%d" % (iBin)):
                self.modelBuilder.out.var("K1Bin%d" % (iBin)).setRange(0.0, 1.0/fracSM4e)
                self.modelBuilder.out.var("K1Bin%d" % (iBin)).setConstant(False)
            else :
                self.modelBuilder.doVar("K1Bin%d[1.0,%s,%s]" % (iBin, 0.0, 1.0/fracSM4e))

            if self.modelBuilder.out.var("K2Bin%d" % (iBin)):
                self.modelBuilder.out.var("K2Bin%d" % (iBin)).setRange(0.0, (1.0-fracSM4e)/fracSM4mu)
                self.modelBuilder.out.var("K2Bin%d" % (iBin)).setConstant(False)
            else :
                self.modelBuilder.doVar("K2Bin%d[1.0,%s,%s]" % (iBin, 0.0, (1.0-fracSM4e)/fracSM4mu))

            if iBin>=0:
                POIs+="SigmaBin%d,"%iBin
                POIs+="K1Bin%d,"%iBin
                POIs+="K2Bin%d,"%iBin
                if self.debug>0:print "Added Bin%d to the POIs"%iBin

        poiNames=[]
        if self.modelBuilder.out.var("MH"):
            if len(self.MHRange) == 2:
                print 'MH will be left floating within', self.MHRange[0], 'and', self.MHRange[1]
                self.modelBuilder.out.var("MH").setRange(float(self.MHRange[0]),float(self.MHRange[1]))
                self.modelBuilder.out.var("MH").setConstant(False)
                poiNames += [ 'MH' ]
            else:
                print 'MH will be assumed to be', self.defaultMH
                self.modelBuilder.out.var("MH").removeRange()
                self.modelBuilder.out.var("MH").setVal(self.defaultMH)
        else:
            if len(self.MHRange) == 2:
                print 'MH will be left floating within', self.MHRange[0], 'and', self.MHRange[1]
                self.modelBuilder.doVar("MH[%s,%s]" % (self.MHRange[0],self.MHRange[1]))
                poiNames += [ 'MH' ]
            else:
                print 'MH (not there before) will be assumed to be', self.defaultMH
                self.modelBuilder.doVar("MH[%g]" % self.defaultMH)
        for poi in poiNames:
            POIs += "%s,"%poi
        POIs = POIs[:-1] # remove last comma
        self.modelBuilder.doSet("POI",POIs)
        self.setup()

    def setup(self):        
        for iBin in range(0,self.nBin):
            self.modelBuilder.factory_('expr::Sigma_trueH4eBin%d("@0*@1*@2", SigmaBin%d, fracSM4eBin%d, K1Bin%d)' % (iBin,iBin,iBin,iBin))
            self.modelBuilder.factory_('expr::Sigma_trueH4muBin%d("@0*(1.0-@1*@2)*@3*@4/(1.0-@1)", SigmaBin%d, fracSM4eBin%d, K1Bin%d, K2Bin%d, fracSM4muBin%d)' % (iBin,iBin,iBin,iBin,iBin,iBin) )
            self.modelBuilder.factory_('expr::Sigma_trueH2e2muBin%d("@0*(1.0-@1*@2)*(1.0-@3*@4/(1.0-@1))", SigmaBin%d, fracSM4eBin%d, K1Bin%d, K2Bin%d, fracSM4muBin%d)' % (iBin,iBin,iBin,iBin,iBin,iBin) )
            self.modelBuilder.factory_('expr::Sigma_trueZ4eBin%d("@0*@1*@2", SigmaBin%d, fracSM4eBin%d, K1Bin%d)' % (iBin,iBin,iBin,iBin))
            self.modelBuilder.factory_('expr::Sigma_trueZ4muBin%d("@0*(1.0-@1*@2)*@3*@4/(1.0-@1)", SigmaBin%d, fracSM4eBin%d, K1Bin%d, K2Bin%d, fracSM4muBin%d)' % (iBin,iBin,iBin,iBin,iBin,iBin) )
            self.modelBuilder.factory_('expr::Sigma_trueZ2e2muBin%d("@0*(1.0-@1*@2)*(1.0-@3*@4/(1.0-@1))", SigmaBin%d, fracSM4eBin%d, K1Bin%d, K2Bin%d, fracSM4muBin%d)' % (iBin,iBin,iBin,iBin,iBin,iBin) )



    def getYieldScale(self,bin,process):
        if not self.DC.isSignal[process]: return 1
        Processes = []
        for Boson in ['H', 'Z']:
            for iBin in range(0,self.nBin):
                for channel in ['4e', '4mu', '2e2mu']:       
                    Processes += ['true'+Boson+channel+'Bin'+str(iBin)]
        if process in Processes: return 'Sigma_'+process
        else: return 1

class DifferentialFiducialK( PhysicsModel ):
    ''' Model used to unfold differential distributions for Fiducial cross-section model for both H->4l and Z->4l'''

    def __init__(self):
        PhysicsModel.__init__(self)
        self.muBinRange=[0.,5.]
        self.nBin=4
        self.MHRange=[20.0,200.0]
        self.defautMH=125.0
        self.fState="inclusive"
        self.debug=1

    def setPhysicsOptions(self,physOptions):
        if self.debug>0:print "Setting PhysicsModel Options"
        for po in physOptions:
            if po.startswith("range="):
                self.muBinRange=po.replace("range=","").split(",")
                if len(self.muBinRange)!=2:
                    raise RunTimeError, "muBinRange require minimal and maximal values: range=min,max"
                if self.debug>0:print "New muBinRange is ", self.muBinRange
            if po.startswith("higgsMassRange="):
                if self.debug>0: print "setting MHRange floating:",po.replace("higgsMassRange=","").split(",")
                self.MHRange=po.replace("higgsMassRange=","").split(",")
                #checks
                if len(self.MHRange) != 2:
                    raise RuntimeError, "MHRange definition requires two extrema: higgsMassRange=min,max"
                elif float(self.MHRange[0]) >= float(self.MHRange[1]):
                    raise RuntimeError, "Extrema for MH defined with inverterd order. Second must be larger the first"
            if po.startswith("mass="):
                self.defaultMH=float( po.replace('mass=','') )
            if po.startswith("nBin="):
                self.nBin=int(po.replace("nBin=",""))
                if self.debug>0:print "new n. of bins is ",self.nBin
            if po.startswith("fState="):
                self.fState=po.replace("fState=","")
                print(self.fState)
                if self.debug>0:print "Final state is ",self.fState
            #verbose
            if po.startswith("verbose"):
                self.debug = 1


    def doParametersOfInterest(self):
        POIs=""
        if self.debug>0:print "Setting pois"

        for iBin in range(0,self.nBin):
            # get values from the workspace
            if(self.fState!='4mu'):
                fracSM4e = self.modelBuilder.out.var("fracSM4eBin%d" % (iBin)).getVal()
            if(self.fState!='4e'):
                fracSM4mu = self.modelBuilder.out.var("fracSM4muBin%d" % (iBin)).getVal()
            SigmaBin = self.modelBuilder.out.var("SigmaBin%d" % (iBin)).getVal()

            if self.modelBuilder.out.var("muBin%d" % (iBin)):
                self.modelBuilder.out.var("muBin%d" % (iBin)).setRange(self.muBinRange[0], self.muBinRange[1])
                self.modelBuilder.out.var("muBin%d" % (iBin)).setConstant(False)
            else :
                self.modelBuilder.doVar("muBin%d[1, %s,%s]" % (iBi, self.muBinRange[0],self.muBinRange[1]))

            if iBin>=0:
                POIs+="muBin%d,"%iBin
                print(self.fState)
                if self.debug>0:print "Added Bin%d to the POIs"%iBin

        poiNames=[]
        if self.modelBuilder.out.var("MH"):
            if len(self.MHRange) == 2:
                print 'MH will be left floating within', self.MHRange[0], 'and', self.MHRange[1]
                self.modelBuilder.out.var("MH").setRange(float(self.MHRange[0]),float(self.MHRange[1]))
                self.modelBuilder.out.var("MH").setConstant(False)
                poiNames += [ 'MH' ]
            else:
                print 'MH will be assumed to be', self.defaultMH
                self.modelBuilder.out.var("MH").removeRange()
                self.modelBuilder.out.var("MH").setVal(self.defaultMH)
        else:
            if len(self.MHRange) == 2:
                print 'MH will be left floating within', self.MHRange[0], 'and', self.MHRange[1]
                self.modelBuilder.doVar("MH[%s,%s]" % (self.MHRange[0],self.MHRange[1]))
                poiNames += [ 'MH' ]
            else:
                print 'MH (not there before) will be assumed to be', self.defaultMH
                self.modelBuilder.doVar("MH[%g]" % self.defaultMH)
        for poi in poiNames:
            POIs += "%s,"%poi
        POIs = POIs[:-1] # remove last comma
        self.modelBuilder.doSet("POI",POIs)
        self.setup()

    def setup(self):        
        for iBin in range(0,self.nBin):
            if (self.fState=='4e'):
                self.modelBuilder.factory_('expr::Sigma_trueH4eBin%d("@0*@1*@2", muBin%d, SigmaBin%d, fracSM4eBin%d)' % (iBin,iBin,iBin,iBin))
            elif (self.fState=='4mu'):
                self.modelBuilder.factory_('expr::Sigma_trueH4muBin%d("@0*@1*@2", muBin%d, SigmaBin%d, fracSM4muBin%d)' % (iBin,iBin,iBin,iBin))
            elif (self.fState=='2e2mu'):
                self.modelBuilder.factory_('expr::Sigma_trueH2e2muBin%d("@0*@1*(1.0-@2-@3)", muBin%d, SigmaBin%d, fracSM4eBin%d, fracSM4muBin%d)' % (iBin,iBin,iBin,iBin,iBin) )
            else:
                self.modelBuilder.factory_('expr::Sigma_trueH4eBin%d("@0*@1*@2", muBin%d, SigmaBin%d, fracSM4eBin%d)' % (iBin,iBin,iBin,iBin))
                self.modelBuilder.factory_('expr::Sigma_trueH4muBin%d("@0*@1*@2", muBin%d, SigmaBin%d, fracSM4muBin%d)' % (iBin,iBin,iBin,iBin))
                self.modelBuilder.factory_('expr::Sigma_trueH2e2muBin%d("@0*@1*(1.0-@2-@3)", muBin%d, SigmaBin%d, fracSM4eBin%d, fracSM4muBin%d)' % (iBin,iBin,iBin,iBin,iBin) )

    def getYieldScale(self,bin,process):
        if not self.DC.isSignal[process]: return 1
        if (self.fState=='4e'):
            fStates = ['4e']
        elif (self.fState=='4mu'):
            fStates = ['4mu']
        elif (self.fState=='2e2mu'):
            fStates = ['2e2mu']
        else:
            fStates = ['4e', '4mu', '2e2mu']

        Processes = []
        Boson = 'H'
        for iBin in range(0,self.nBin):
            for channel in fStates:       
                Processes += ['true'+Boson+channel+'Bin'+str(iBin)]
        if process in Processes: return 'Sigma_'+process
        else: return 1


class H4lZ4lInclusiveFiducialRatio( PhysicsModel ):
    ''' Model used to unfold differential distributions '''

    def __init__(self):
        PhysicsModel.__init__(self)
        self.SigmaRange=[0.,10]
        self.RatioSigmaRange=[0.,10]
        self.MHRange=[115.,130.]
        self.DeltaMHmZRange=[30.,40.]
        self.defaultMH = 125.0
        self.defaultDeltaMHmZ = 33.8124
        self.fixMH = False
        self.fixDeltaMHmZ = False
        self.debug=1

    def setPhysicsOptions(self,physOptions):
        if self.debug>0:print "Setting PhysicsModel Options"
        for po in physOptions:
            if po.startswith("SigmaRange="):
                self.SigmaRange=po.replace("SigmaRange=","").split(":")
                if len(self.SigmaRange)!=2:
                    raise RunTimeError, "SigmaRange require minimal and maximal values: SigmaRange=min:max"
                if self.debug>0:print "new SigmaRange is ", self.SigmaRange
            if po.startswith("RatioSigmaRange="):
                self.RatioSigmaRange=po.replace("RatioSigmaRange=","").split(":")
                if len(self.RatioSigmaRange)!=2:
                    raise RunTimeError, "RatioSigmaRange require minimal and maximal values: RatioSigmaRange=min:max"
                if self.debug>0:print "new RatioSigmaRange is ", self.RatioSigmaRange
            if po.startswith("MHRange="):
                if self.debug>0: print "setting MH mass range floating:",po.replace("MHRange=","").split(":")
                self.MHRange=po.replace("MHRange=","").split(",")
                #checks
                if len(self.MHRange) != 2:
                    raise RuntimeError, "MH range definition requires two extrema"
                elif float(self.MHRange[0]) >= float(self.MHRange[1]):
                    raise RuntimeError, "Extrema for MH mass range defined with inverterd order. Second must be larger the first"
            if po.startswith("DeltaMHmZRange="):
                if self.debug>0: print "setting MH-MZ mass range floating:",po.replace("DeltaMHmZRange=","").split(":")
                self.DeltaMHmZRange=po.replace("DeltaMHmZRange=","").split(",")
                #checks
                if len(self.DeltaMHmZRange) != 2:
                    raise RuntimeError, "MH-MZ range definition requires two extrema"
                elif float(self.DeltaMHmZRange[0]) >= float(self.DeltaMHmZRange[1]):
                    raise RuntimeError, "Extrema for MH-MZ mass range defined with inverterd order. Second must be larger the first"
            if po.startswith("defaultMH="):
                self.defaultMH=float( po.replace('defaultMH=','') )
            if po.startswith("defaultDeltaMHmZ="):
                self.defaultDeltaMHmZ=float( po.replace('defaultDeltaMHmZ=','') )
            if po.startswith("fixMH"):
                self.fixMH = True
            if po.startswith("fixDeltaMHmZ"):
                self.fixDeltaMHmZ = True
            #verbose
            if po.startswith("verbose"):
                self.debug = 1

    def doParametersOfInterest(self):
        POIs=""
        if self.debug>0:print "Setting pois"

        if self.modelBuilder.out.var("SigmaH"):
            self.modelBuilder.out.var("SigmaH").setRange(self.SigmaRange[0], self.SigmaRange[1])
        else: 
            self.modelBuilder.doVar("SigmaH[1,%s,%s]" % (self.SigmaRange[0], self.SigmaRange[1]))
        if self.modelBuilder.out.var("SigmaH4e"):
            self.modelBuilder.out.var("SigmaH4e").setRange(self.SigmaRange[0], self.SigmaRange[1])
        else:        
            self.modelBuilder.doVar("SigmaH4e[1,%s,%s]" % (self.SigmaRange[0], self.SigmaRange[1]))
        if self.modelBuilder.out.var("SigmaH4mu"):
            self.modelBuilder.out.var("SigmaH4mu").setRange(self.SigmaRange[0], self.SigmaRange[1])
        else:
            self.modelBuilder.doVar("SigmaH4mu[1,%s,%s]" % (self.SigmaRange[0], self.SigmaRange[1]))
        if self.modelBuilder.out.var("RatioSigmaHoZ"):
            self.modelBuilder.out.var("RatioSigmaHoZ").setRange(self.RatioSigmaRange[0], self.RatioSigmaRange[1])
        else:
            self.modelBuilder.doVar("RatioSigmaHoZ[1,%s,%s]" % (self.RatioSigmaRange[0], self.RatioSigmaRange[1]))
        if self.modelBuilder.out.var("RatioSigmaHoZ4e"):
            self.modelBuilder.out.var("RatioSigmaHoZ4e").setRange(self.RatioSigmaRange[0], self.RatioSigmaRange[1])
        else:
            self.modelBuilder.doVar("RatioSigmaHoZ4e[1,%s,%s]" % (self.RatioSigmaRange[0], self.RatioSigmaRange[1]))
        if self.modelBuilder.out.var("RatioSigmaHoZ4mu"):
            self.modelBuilder.out.var("RatioSigmaHoZ4mu").setRange(self.RatioSigmaRange[0], self.RatioSigmaRange[1])
        else:
            self.modelBuilder.doVar("RatioSigmaHoZ4mu[1,%s,%s]" % (self.RatioSigmaRange[0], self.RatioSigmaRange[1]))

        POIs+="SigmaH,"
        POIs+="SigmaH4e,"
        POIs+="SigmaH4mu,"
        POIs+="RatioSigmaHoZ,"
        POIs+="RatioSigmaHoZ4e,"
        POIs+="RatioSigmaHoZ4mu,"

        # --- Other parameters ----
        poiNames=[]
        # set Parameter MH
        if self.modelBuilder.out.var("MH"):
            if len(self.MHRange) == 2:
                print 'MH will be left floating within', self.MHRange[0], 'and', self.MHRange[1]
                self.modelBuilder.out.var("MH").setRange(float(self.MHRange[0]),float(self.MHRange[1]))
                self.modelBuilder.out.var("MH").setVal(self.defaultMH)
                self.modelBuilder.out.var("MH").setConstant(False)
                poiNames += [ 'MH' ]
            else:
                print 'MH will be assumed to be', self.defaultMH
                self.modelBuilder.out.var("MH").removeRange()
                self.modelBuilder.out.var("MH").setVal(self.defaultMH)
                self.modelBuilder.out.var("MH").setConstant(True)
        else:
            if len(self.MHRange) == 2:
                print 'MH will be left floating within', self.MHRange[0], 'and', self.MHRange[1]
                self.modelBuilder.doVar("MH[%s,%s]" % (self.MHRange[0],self.MHRange[1]))
                self.modelBuilder.out.var("MH").setVal(self.defaultMH)
                poiNames += [ 'MH' ]
            else:
                print 'MH (not there before) will be assumed to be', self.defaultMH
                self.modelBuilder.doVar("MH[%g]" % self.defaultMH)
        # set Parameter DeltaMHmZ 
        if self.modelBuilder.out.var("DeltaMHmZ"):
            if len(self.DeltaMHmZRange) == 2:
                print 'DeltaMHmZ will be left floating within', self.DeltaMHmZRange[0], 'and', self.DeltaMHmZRange[1]
                self.modelBuilder.out.var("DeltaMHmZ").setRange(float(self.DeltaMHmZRange[0]),float(self.DeltaMHmZRange[1]))
                self.modelBuilder.out.var("DeltaMHmZ").setVal(self.defaultDeltaMHmZ)
                self.modelBuilder.out.var("DeltaMHmZ").setConstant(False)
                poiNames += [ 'DeltaMHmZ' ]
            else:
                print 'DeltaMHmZ will be assumed to be', self.defaultDeltaMHmZ
                self.modelBuilder.out.var("DeltaMHmZ").removeRange()
                self.modelBuilder.out.var("DeltaMHmZ").setVal(self.defaultDeltaMHmZ)
                self.modelBuilder.out.var("DeltaMHmZ").setConstant(True)
        else:
            if len(self.DeltaMHmZRange) == 2:
                print 'DeltaMHmZ will be left floating within', self.DeltaMHmZRange[0], 'and', self.DeltaMHmZRange[1]
                self.modelBuilder.doVar("DeltaMHmZ[%s,%s]" % (self.DeltaMHmZRange[0],self.DeltaMHmZRange[1]))
                self.modelBuilder.out.var("DeltaMHmZ").setVal(self.defaultDeltaMHmZ)
                poiNames += [ 'DeltaMHmZ' ]
            else:
                print 'DeltaMHmZ (not there before) will be assumed to be', self.defaultDeltaMHmZ
                self.modelBuilder.doVar("DeltaMHmZ[%g]" % self.defaultDeltaMHmZ)

        if (self.fixMH):
            self.modelBuilder.out.var("MH").setVal(self.defaultMH)
            self.modelBuilder.out.var("MH").setConstant(True)
            for ixx in range(poiNames.count('MH')): poiNames.remove('MH')
        if (self.fixDeltaMHmZ):
            self.modelBuilder.out.var("DeltaMHmZ").setVal(self.defaultDeltaMHmZ)
            self.modelBuilder.out.var("DeltaMHmZ").setConstant(True)
            for ixx in range(poiNames.count('DeltaMHmZ')): poiNames.remove('DeltaMHmZ')

        for poi in poiNames:
            POIs += "%s,"%poi
        POIs = POIs[:-1] # remove last comma
        self.modelBuilder.doSet("POI",POIs)
        print "set up pois"
        self.setup()

    def setup(self):
        self.modelBuilder.factory_('expr::sigma_trueH4e("@0", SigmaH4e)')
        self.modelBuilder.factory_('expr::sigma_trueH4mu("@0", SigmaH4mu)')
        self.modelBuilder.factory_('expr::sigma_trueH2e2mu("(@0-@1-@2)", SigmaH, SigmaH4e, SigmaH4mu)')
        self.modelBuilder.factory_('expr::sigma_trueZ4e("(@0/@1)", SigmaH4e, RatioSigmaHoZ4e)')
        self.modelBuilder.factory_('expr::sigma_trueZ4mu("(@0/@1)", SigmaH4mu, RatioSigmaHoZ4mu)')
        self.modelBuilder.factory_('expr::sigma_trueZ2e2mu("(@0/@1-@2/@3-@4/@5)", SigmaH, RatioSigmaHoZ, SigmaH4e, RatioSigmaHoZ4e, SigmaH4mu, RatioSigmaHoZ4mu)')

    def getYieldScale(self,bin,process):
        if not self.DC.isSignal[process]: return 1
        if process in [ "trueH4e", "trueH4mu","trueH2e2mu","trueZ4e", "trueZ4mu","trueZ2e2mu"]:
            return "sigma_"+process
        else: 
            return 1


class H4lZ4lInclusiveFiducialRatioV2( PhysicsModel ):
    ''' Model used to unfold differential distributions '''

    def __init__(self):
        PhysicsModel.__init__(self)
        self.SigmaRange=[0.,10]
        self.RatioSigmaRange=[0.,10]
        self.MHRange=[115.,130.]
        self.DeltaMHmZRange=[30.,40.]
        self.defaultMH = 125.0
        self.defaultDeltaMHmZ = 33.8124
        self.fixMH = False
        self.fixDeltaMHmZ = False
        self.debug=1

    def setPhysicsOptions(self,physOptions):
        if self.debug>0:print "Setting PhysicsModel Options"
        for po in physOptions:
            if po.startswith("SigmaRange="):
                self.SigmaRange=po.replace("SigmaRange=","").split(":")
                if len(self.SigmaRange)!=2:
                    raise RunTimeError, "SigmaRange require minimal and maximal values: SigmaRange=min:max"
                if self.debug>0:print "new SigmaRange is ", self.SigmaRange
            if po.startswith("RatioSigmaRange="):
                self.RatioSigmaRange=po.replace("RatioSigmaRange=","").split(":")
                if len(self.RatioSigmaRange)!=2:
                    raise RunTimeError, "RatioSigmaRange require minimal and maximal values: RatioSigmaRange=min:max"
                if self.debug>0:print "new RatioSigmaRange is ", self.RatioSigmaRange
            if po.startswith("MHRange="):
                if self.debug>0: print "setting MH mass range floating:",po.replace("MHRange=","").split(":")
                self.MHRange=po.replace("MHRange=","").split(",")
                #checks
                if len(self.MHRange) != 2:
                    raise RuntimeError, "MH range definition requires two extrema"
                elif float(self.MHRange[0]) >= float(self.MHRange[1]):
                    raise RuntimeError, "Extrema for MH mass range defined with inverterd order. Second must be larger the first"
            if po.startswith("DeltaMHmZRange="):
                if self.debug>0: print "setting MH-MZ mass range floating:",po.replace("DeltaMHmZRange=","").split(":")
                self.DeltaMHmZRange=po.replace("DeltaMHmZRange=","").split(",")
                #checks
                if len(self.DeltaMHmZRange) != 2:
                    raise RuntimeError, "MH-MZ range definition requires two extrema"
                elif float(self.DeltaMHmZRange[0]) >= float(self.DeltaMHmZRange[1]):
                    raise RuntimeError, "Extrema for MH-MZ mass range defined with inverterd order. Second must be larger the first"
            if po.startswith("defaultMH="):
                self.defaultMH=float( po.replace('defaultMH=','') )
            if po.startswith("defaultDeltaMHmZ="):
                self.defaultDeltaMHmZ=float( po.replace('defaultDeltaMHmZ=','') )
            if po.startswith("fixMH"):
                self.fixMH = True
            if po.startswith("fixDeltaMHmZ"):
                self.fixDeltaMHmZ = True
            #verbose
            if po.startswith("verbose"):
                self.debug = 1

    def doParametersOfInterest(self):
        POIs=""
        if self.debug>0:print "Setting pois"
        if self.modelBuilder.out.var("SigmaH4e"):
            self.modelBuilder.out.var("SigmaH4e").setRange(self.SigmaRange[0], self.SigmaRange[1])
        else:
            self.modelBuilder.doVar("SigmaH4e[1,%s,%s]" % (self.SigmaRange[0], self.SigmaRange[1]))
        if self.modelBuilder.out.var("SigmaH4mu"):
            self.modelBuilder.out.var("SigmaH4mu").setRange(self.SigmaRange[0], self.SigmaRange[1])
        else:
            self.modelBuilder.doVar("SigmaH4mu[1,%s,%s]" % (self.SigmaRange[0], self.SigmaRange[1]))
        if self.modelBuilder.out.var("SigmaH2e2mu"):
            self.modelBuilder.out.var("SigmaH2e2mu").setRange(self.SigmaRange[0], self.SigmaRange[1])
        else:
            self.modelBuilder.doVar("SigmaH2e2mu[1,%s,%s]" % (self.SigmaRange[0], self.SigmaRange[1]))
        if self.modelBuilder.out.var("RatioSigmaH4e"):
            self.modelBuilder.out.var("RatioSigmaH4e").setRange(self.RatioSigmaRange[0], self.RatioSigmaRange[1])
        else:
            self.modelBuilder.doVar("RatioSigmaHoZ4e[1,%s,%s]" % (self.RatioSigmaRange[0], self.RatioSigmaRange[1]))
        if self.modelBuilder.out.var("RatioSigmaH4mu"):
            self.modelBuilder.out.var("RatioSigmaH4mu").setRange(self.RatioSigmaRange[0], self.RatioSigmaRange[1])
        else:
            self.modelBuilder.doVar("RatioSigmaHoZ4mu[1,%s,%s]" % (self.RatioSigmaRange[0], self.RatioSigmaRange[1]))
        if self.modelBuilder.out.var("RatioSigmaH2e2mu"):
            self.modelBuilder.out.var("RatioSigmaH2e2mu").setRange(self.RatioSigmaRange[0], self.RatioSigmaRange[1])
        else:
            self.modelBuilder.doVar("RatioSigmaHoZ2e2mu[1,%s,%s]" % (self.RatioSigmaRange[0], self.RatioSigmaRange[1]))

        POIs+="SigmaH4e,"
        POIs+="SigmaH4mu,"
        POIs+="SigmaH2e2mu,"
        POIs+="RatioSigmaHoZ4e,"
        POIs+="RatioSigmaHoZ4mu,"
        POIs+="RatioSigmaHoZ2e2mu,"

        # --- Other parameters ----
        poiNames=[]
        # set Parameter MH
        if self.modelBuilder.out.var("MH"):
            if len(self.MHRange) == 2:
                print 'MH will be left floating within', self.MHRange[0], 'and', self.MHRange[1]
                self.modelBuilder.out.var("MH").setRange(float(self.MHRange[0]),float(self.MHRange[1]))
                self.modelBuilder.out.var("MH").setVal(self.defaultMH)
                self.modelBuilder.out.var("MH").setConstant(False)
                poiNames += [ 'MH' ]
            else:
                print 'MH will be assumed to be', self.defaultMH
                self.modelBuilder.out.var("MH").removeRange()
                self.modelBuilder.out.var("MH").setVal(self.defaultMH)
                self.modelBuilder.out.var("MH").setConstant(True)
        else:
            if len(self.MHRange) == 2:
                print 'MH will be left floating within', self.MHRange[0], 'and', self.MHRange[1]
                self.modelBuilder.doVar("MH[%s,%s]" % (self.MHRange[0],self.MHRange[1]))
                self.modelBuilder.out.var("MH").setVal(self.defaultMH)
                poiNames += [ 'MH' ]
            else:
                print 'MH (not there before) will be assumed to be', self.defaultMH
                self.modelBuilder.doVar("MH[%g]" % self.defaultMH)
        # set Parameter DeltaMHmZ
        if self.modelBuilder.out.var("DeltaMHmZ"):
            if len(self.DeltaMHmZRange) == 2:
                print 'DeltaMHmZ will be left floating within', self.DeltaMHmZRange[0], 'and', self.DeltaMHmZRange[1]
                self.modelBuilder.out.var("DeltaMHmZ").setRange(float(self.DeltaMHmZRange[0]),float(self.DeltaMHmZRange[1]))
                self.modelBuilder.out.var("DeltaMHmZ").setVal(self.defaultDeltaMHmZ)
                self.modelBuilder.out.var("DeltaMHmZ").setConstant(False)
                poiNames += [ 'DeltaMHmZ' ]
            else:
                print 'DeltaMHmZ will be assumed to be', self.defaultDeltaMHmZ
                self.modelBuilder.out.var("DeltaMHmZ").removeRange()
                self.modelBuilder.out.var("DeltaMHmZ").setVal(self.defaultDeltaMHmZ)
                self.modelBuilder.out.var("DeltaMHmZ").setConstant(True)
        else:
            if len(self.DeltaMHmZRange) == 2:
                print 'DeltaMHmZ will be left floating within', self.DeltaMHmZRange[0], 'and', self.DeltaMHmZRange[1]
                self.modelBuilder.doVar("DeltaMHmZ[%s,%s]" % (self.DeltaMHmZRange[0],self.DeltaMHmZRange[1]))
                self.modelBuilder.out.var("DeltaMHmZ").setVal(self.defaultDeltaMHmZ)
                poiNames += [ 'DeltaMHmZ' ]
            else:
                print 'DeltaMHmZ (not there before) will be assumed to be', self.defaultDeltaMHmZ
                self.modelBuilder.doVar("DeltaMHmZ[%g]" % self.defaultDeltaMHmZ)

        if (self.fixMH):
            self.modelBuilder.out.var("MH").setVal(self.defaultMH)
            self.modelBuilder.out.var("MH").setConstant(True)
            for ixx in range(poiNames.count('MH')): poiNames.remove('MH')
        if (self.fixDeltaMHmZ):
            self.modelBuilder.out.var("DeltaMHmZ").setVal(self.defaultDeltaMHmZ)
            self.modelBuilder.out.var("DeltaMHmZ").setConstant(True)
            for ixx in range(poiNames.count('DeltaMHmZ')): poiNames.remove('DeltaMHmZ')

        for poi in poiNames:
            POIs += "%s,"%poi
        POIs = POIs[:-1] # remove last comma
        self.modelBuilder.doSet("POI",POIs)
        print "set up pois"
        self.setup()

    def setup(self):
        self.modelBuilder.factory_('expr::sigma_trueH4e("@0", SigmaH4e)')
        self.modelBuilder.factory_('expr::sigma_trueH4mu("@0", SigmaH4mu)')
        self.modelBuilder.factory_('expr::sigma_trueH2e2mu("@0", SigmaH2e2mu)')
        self.modelBuilder.factory_('expr::sigma_trueZ4e("(@0/@1)", SigmaH4e, RatioSigmaHoZ4e)')
        self.modelBuilder.factory_('expr::sigma_trueZ4mu("(@0/@1)", SigmaH4mu, RatioSigmaHoZ4mu)')
        self.modelBuilder.factory_('expr::sigma_trueZ2e2mu("(@0/@1)", SigmaH2e2mu, RatioSigmaHoZ2e2mu)')

    def getYieldScale(self,bin,process):
        if not self.DC.isSignal[process]: return 1
        if process in [ "trueH4e", "trueH4mu","trueH2e2mu","trueZ4e", "trueZ4mu","trueZ2e2mu"]:
            return "sigma_"+process
        else:
            return 1

class TrilinearHiggs(PhysicsModel):
    "Float independently cross sections and branching ratios"
    # def __init__(self):
    #     PhysicsModel.__init__(self)
    #     # SMLikeHiggsModel.__init__(self) # not using 'super(x,self).__init__' since I don't understand it
    #     self.mHRange = []
    #     self.poiNames = []

    def __init__(self):
        PhysicsModel.__init__(self)
        self.nBin=4
        self.MHRange=[20.0,200.0]
        self.defautMH=125.38

    def setPhysicsOptions(self,physOptions):
        for po in physOptions:
            if po.startswith("higgsMassRange="):
                self.mHRange = po.replace("higgsMassRange=","").split(",")
                if len(self.mHRange) != 2:
                    raise RuntimeError, "Higgs mass range definition requires two extrema"
                elif float(self.mHRange[0]) >= float(self.mHRange[1]):
                    raise RuntimeError, "Extrema for Higgs mass range defined with inverterd order. Second must be larger the first"

    def doParametersOfInterest(self):
        """Create POI and other parameters, and define the POI set."""
        POIs=""
        # trilinear Higgs couplings modified 
        self.modelBuilder.doVar("k_lambda[1,-20.,20.]")
        # self.poiNames="k_lambda"
        POIs = "k_lambda"
        self.modelBuilder.doSet("POI",POIs)
        print POIs
        self.setup()

    def setup(self):
        # Let's start with ggH
        #Use inclusive value for ggH: EWK reweighting tool not available. Taken directly from arXiv:1607.04251
        proc = "ggH"
        C1_ggH = 0.0066
        C1_map = {}
        for i in range(5): # equivalent to nBins above
            C1_map["ggH_gen%g"%i] = C1_ggH
 
        #Define dZH constant variable
        dZH = -1.536e-3

        #Loop over processes*gen bins in map to define how cross-section scales
        # for proc in C1_map:
        #   self.modelBuilder.factory_("expr::XSscal_%s(\"(1+@0*%g+%g)/((1-(@0*@0-1)*%g)*(1+%g+%g))\",k_lambda)"%(proc,C1_map[proc],dZH,dZH,C1_map[proc],dZH))
        # For the moment ggH scaling only
        self.modelBuilder.factory_("expr::XSscal_%s(\"(1+@0*%g+%g)/((1-(@0*@0-1)*%g)*(1+%g+%g))\",k_lambda)"%(proc,C1_ggH,dZH,dZH,C1_ggH,dZH))

        #Scaling @ decay: define expression for how BR scales as function of klambda: h->gammagamma
        #Use following parameters taken directly from: arXiv:1607.04251
        # C1_hgg = 0.0049
        C1_hzz = 0.0083 # hzz4l
        C1_tot = 2.5e-3
        self.modelBuilder.factory_("expr::BRscal_hzz(\"1+(((@0-1)*(%g-%g))/(1+(@0-1)*%g))\",k_lambda)"%(C1_hzz,C1_tot,C1_tot))
        self.modelBuilder.factory_('expr::XSBRscal_ggH_hzz(\"(@0*@1)\", XSscal_ggH, BRscal_hzz)')
        
        # print self.poiNames
        # self.modelBuilder.doSet("POI",self.poiNames)

    def getYieldScale(self,bin,process):
        if not self.DC.isSignal[process]: return 1
        name = "XSBRscal_ggH_hzz" #% (production,decay)
        # self.modelBuilder.factory_('expr::%s(\"(@0*@1)\", XSscal_ggH, BRscal_hzz)'%(name))
        return name

        # Processes = []
        # Boson = 'H'
        # for iBin in range(5): #,self.nBin):
        #     for channel in fStates:       
        #         Processes += ['ggH_gen'+str(iBin)+'_hzz']
        # if process in Processes: 

        #     return 'Sigma_'+process
        # else: return 1

    # def getHiggsSignalYieldScale(self,production,decay):
        
    #     #XSBR
    #     name = "XSBRscal_%s_%s" % (production,decay)
    #     #Name has not been defined in doParametersOfInterest: combine XS + BR
    #     if self.modelBuilder.out.function(name) == None:
    #       #XS
    #       if self.modelBuilder.out.function( "XSscal_%s"%(production) ) == None:
    #         print "DEBUG: proc not given XS scaling"
    #         raise RuntimeError, "Production mode %s not supported"%production
    #       else:
    #         XSscal = "XSscal_%s_%s"%(production,decay)
    #       #BR
    #       if self.modelBuilder.out.function( "BRscal_%s"%(decay) ) == None:
    #         print "DEBUG: proc not given BR scaling"
    #         raise RuntimeError, "Decay mode %s not supported"%decay
    #       else:
    #         BRscal = "BRscal_%s"%decay
    #       #XSBR
    #       self.modelBuilder.factory_('expr::%s(\"(@0*@1)\", XSscal_%s, BRscal_%s)'%(name,production,decay))
    #       print '[LHC-CMS Trilinear]', name, ": ", self.modelBuilder.out.function(name).Print("")
    #     return name
                     
inclusiveFiducial=InclusiveFiducial()
inclusiveFiducialV2=InclusiveFiducialV2()
inclusiveFiducialV3=InclusiveFiducialV3()

differentialFiducial=DifferentialFiducial()
differentialFiducialV2=DifferentialFiducialV2()
differentialFiducialV3=DifferentialFiducialV3()
differentialFiducialK=DifferentialFiducialK()

h4lZ4lInclusiveFiducialRatio=H4lZ4lInclusiveFiducialRatio()
h4lZ4lInclusiveFiducialRatioV2=H4lZ4lInclusiveFiducialRatioV2()


trilinearHiggs = TrilinearHiggs()

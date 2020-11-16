from HiggsAnalysis.CombinedLimit.PhysicsModel import *

class InclusiveFiducialV2( PhysicsModel ):
        ''' Model used to unfold differential distributions '''

        def __init__(self):
                PhysicsModel.__init__(self)
                self.Range=[0.,10.0]
                self.fracRange=[0.,0.5]                
                self.mHRange=[]
                self.debug=1
                self.mass=0
                #regularization is in the datacard
                #self.deltaReg=0
                #self.regBins=""

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
                
                self.modelBuilder.doVar("r2e2mu[1,%s,%s]" % (self.Range[0], self.Range[1]))
                self.modelBuilder.doVar("r4e[1,%s,%s]" % (self.Range[0], self.Range[1]))
                self.modelBuilder.doVar("r4mu[1,%s,%s]" % (self.Range[0], self.Range[1]))
                
                POIs+="r2e2mu,"
                POIs+="r4e,"
                POIs+="r4mu,"
                                        
                # --- Higgs Mass as other parameter ----
#               if self.options.mass != 0:
#                   if self.modelBuilder.out.var("MH"):
#                     self.modelBuilder.out.var("MH").removeRange()
#                     self.modelBuilder.out.var("MH").setVal(self.options.mass)
#                   else:
#                     self.modelBuilder.doVar("MH[%g]" % self.options.mass);
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

             self.modelBuilder.factory_('expr::r_trueH2e2mu("@0", r2e2mu)')
             self.modelBuilder.factory_('expr::r_trueH4e("@0", r4e)')
             self.modelBuilder.factory_('expr::r_trueH4mu("@0", r4mu)')

        def getYieldScale(self,bin,process):

             if not self.DC.isSignal[process]: return 1

             #print "print process"
             #print process
             name = "r_%s" % process
             
             if process in [ "trueH2e2mu", "trueH4e","trueH4mu"]: 
                return name
             else : return 1

            

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
                #regularization is in the datacard
                #self.deltaReg=0
                #self.regBins=""

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
                        #if iBin==self.nBin-1: #for out-of-acceptance bin 
                                #if self.debug>0:print "   and set constant to the value 1 "
                                #self.modelBuilder.out.var("r_Bin%d"%iBin).removeRange()
                                #self.modelBuilder.out.var("r_Bin%d"%iBin).setVal(1)
                                #self.modelBuilder.out.var("r_Bin%d"%iBin).setConstant(True)
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
                self.setup()

        def setup(self):        
               # self.modelBuilder.factory_('expr::%s("@0*@1", %s, frac_%s)' % (name, fiducial, process))                         
               for iBin in range(0,self.nBin):
                 self.modelBuilder.factory_('expr::r_trueH2e2muBin%d("@0", r2e2muBin%d)' % (iBin,iBin))
                 self.modelBuilder.factory_('expr::r_trueH4eBin%d("@0", r4eBin%d)'% (iBin,iBin))
                 self.modelBuilder.factory_('expr::r_trueH4muBin%d("@0", r4muBin%d)'% (iBin,iBin))

        def getYieldScale(self,bin,process):

             if not self.DC.isSignal[process]: return 1

             #print "print process"
             #print process
             name = "fiducial_%s" % process
             
             self.modelBuilder.factory_('expr::%s("@0", r_%s)' % (name, process))
  
             if process in [ "trueH2e2muBin0","trueH4eBin0","trueH4muBin0","trueH2e2muBin1","trueH4eBin1","trueH4muBin1","trueH2e2muBin2","trueH4eBin2","trueH4muBin2","trueH2e2muBin3","trueH4eBin3","trueH4muBin3","trueH2e2muBin4","trueH4eBin4","trueH4muBin4","trueH2e2muBin6","trueH4eBin6","trueH4muBin6","trueH2e2muBin5","trueH4eBin5","trueH4muBin5"]: 
                return name
             else : return 1
               
class DifferentialFiducialV2_4mu( PhysicsModel ):
        ''' Model used to unfold differential distributions '''

        def __init__(self):
                PhysicsModel.__init__(self)
                self.Range=[0.,4]
                self.nBin=4
                self.debug=1
                #regularization is in the datacard
                #self.deltaReg=0
                #self.regBins=""

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
                for iBin in range(0,self.nBin):
                        self.modelBuilder.doVar("r4muBin%d[1,%s,%s]" % (iBin,self.Range[0], self.Range[1]))
                        if iBin>-1:
                          POIs+="r4muBin%d,"%iBin
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
                self.setup()

        def setup(self):
               for iBin in range(0,self.nBin):
                 self.modelBuilder.factory_('expr::r_trueH4muBin%d("@0", r4muBin%d)'% (iBin,iBin))

        def getYieldScale(self,bin,process):

             if not self.DC.isSignal[process]: return 1

             #print "print process"
             #print process
             name = "fiducial_%s" % process

             self.modelBuilder.factory_('expr::%s("@0", r_%s)' % (name, process))

             if process in [ "trueH4muBin0","trueH4muBin1","trueH4muBin2","trueH4muBin3"]:
                return name
             else : return 1


class DifferentialFiducialV2_4e( PhysicsModel ):

        ''' Model used to unfold differential distributions '''

        def __init__(self):
                PhysicsModel.__init__(self)
                self.Range=[0.,4]
                self.nBin=4
                self.debug=1

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
                for iBin in range(0,self.nBin):
                        self.modelBuilder.doVar("r4eBin%d[1,%s,%s]" % (iBin,self.Range[0], self.Range[1]))
                        if iBin>-1:
                          POIs+="r4eBin%d,"%iBin
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
                        POIs += ",%s"%poi
                self.modelBuilder.doSet("POI",POIs)
                self.setup()

        def setup(self):
               for iBin in range(0,self.nBin):
                 self.modelBuilder.factory_('expr::r_trueH4eBin%d("@0", r4eBin%d)'% (iBin,iBin))


        def getYieldScale(self,bin,process):

             if not self.DC.isSignal[process]: return 1

             #print "print process"
             #print process
             name = "fiducial_%s" % process

             self.modelBuilder.factory_('expr::%s("@0", r_%s)' % (name, process))

             if process in [ "trueH4eBin0","trueH4eBin1","trueH4eBin2","trueH4eBin3"]:
                return name
             else : return 1

class DifferentialFiducialV2_2e2mu( PhysicsModel ):
        ''' Model used to unfold differential distributions '''

        def __init__(self):
                PhysicsModel.__init__(self)
                self.Range=[0.,4]
                self.nBin=4
                self.debug=1

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
                for iBin in range(0,self.nBin):
                        self.modelBuilder.doVar("r2e2muBin%d[1,%s,%s]" % (iBin,self.Range[0], self.Range[1]))
                        if iBin>-1:
                          POIs+="r2e2muBin%d,"%iBin
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
                        POIs += ",%s"%poi
                self.modelBuilder.doSet("POI",POIs)
                self.setup()

        def setup(self):
               for iBin in range(0,self.nBin):
                 self.modelBuilder.factory_('expr::r_trueH2e2muBin%d("@0", r2e2muBin%d)'% (iBin,iBin))


        def getYieldScale(self,bin,process):

             if not self.DC.isSignal[process]: return 1

             #print "print process"
             #print process
             name = "fiducial_%s" % process

             self.modelBuilder.factory_('expr::%s("@0", r_%s)' % (name, process))

             if process in [ "trueH2e2muBin0","trueH2e2muBin1","trueH2e2muBin2","trueH2e2muBin3"]:
                return name
             else : return 1


inclusiveFiducialV2=InclusiveFiducialV2()

differentialFiducialV2=DifferentialFiducialV2()

differentialFiducialV2_2e2mu=DifferentialFiducialV2_2e2mu()

differentialFiducialV2_4e=DifferentialFiducialV2_4e()

differentialFiducialV2_4mu=DifferentialFiducialV2_4mu()



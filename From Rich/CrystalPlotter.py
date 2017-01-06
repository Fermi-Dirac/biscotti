from mayavi import mlab
import numpy as np

class Crystal():
    __fname = ''
    __cp  = []   #cell positions
    __ap  = []   #atom positions
    __Ns  = 0    #Atomic Species Count
    __SP  = []   #Atomic Species Names
    __Runs= 0    #Number of Runs

    
    def __init__(self,fname):
        self.loadFile(fname)
        
    def getFname(self):
        return self.__fname
    
    def loadFile(self,fname):
        self.__fname=fname
        with open(fname,'r') as f:
            s = f.readlines()
        sp = [i.split() for i in s]
        cind = [i for i,l in enumerate(s) if 'CELL_PARAMETERS' in l]
        aind = [i for i,l in enumerate(s) if 'ATOMIC_POSITIONS' in l]

        def getNS(ai):
            count = 0
            while True:
                if len(s[ai+count+1].split())==0:
                    break
                count += 1
                return count
        N = getNS(aind[0])
        self.__Ns = N
        self.__SP = np.array(sp[aind[0]+1:aind[0]+1+N])[:,0]
        self.__cp = np.array([np.array(sp[i+1:j-1]) for i,j in zip(cind,aind)]).astype(np.float)
        self.__ap = np.array([np.array(sp[i+1:i+2+N]) for i in aind])[:,:,1:].astype(np.float)
        self.__Runs=len(cind)
    
    def mkCrystalPoints(self,Run,N):
        mult = lambda A,r: (r*A.transpose()).transpose()
        R = lambda h,k,l,A: h*A[0]+k*A[1]+l*A[2]
        Rs = lambda h,k,l,A,r: np.sum(mult(np.array(A),np.array([h,k,l]))+mult(np.array(A),np.array(r)),axis=0)

        Points = lambda CP,SP,N: np.array([Rs(i,j,k,CP,SP) 
                                         for i in range(N) 
                                         for j in range(N) 
                                         for k in range(N)]).transpose()
        count = 1
        
        x,y,z,c=[],[],[],[]
        #x,y,z = Points(self.__cp[Run],np.zeros(3),N)
        #c = np.ones(len(x))*0
        for i in self.__ap[Run]:
            Ix,Iy,Iz = Points(self.__cp[Run],i,N)
            Ic = np.ones(len(Ix))*count
            x = np.concatenate((x,Ix),axis=0)
            y = np.concatenate((y,Iy),axis=0)
            z = np.concatenate((z,Iz),axis=0)
            c = np.concatenate((c,Ic),axis=0)
            count +=1
        #pts = mlab.points3d(x,y,z,c,colormap='jet',scale_factor=1)
        sV = np.ones(len(x))
        return x,y,z,c,sV
    
    def plotCrystal(self,Run,N,scale=.5):
        x,y,z,c,sV = self.mkCrystalPoints(Run,N)
        pts = mlab.quiver3d(x, y, z, sV, sV, sV, scalars=c, mode="sphere"
                            , scale_factor=scale,vmin=0,colormap='jet')
        pts.glyph.color_mode = "color_by_scalar"
        pts.glyph.glyph_source.glyph_source.center = [0,0,0]
        print type(pts)
        mlab.show()
    
    def getNumberRuns(self):
        return self.__Runs
    
    def getPositions(self):
        return self.__cp,self.__ap

#####################
######Interact#######
#####################


from mayavi import mlab
import numpy as np

from traits.api import HasTraits,Int, Range, Instance, Button, on_trait_change
from traitsui.api import View, Item, Group

from mayavi.core.api import PipelineBase
from mayavi.core.ui.api import MayaviScene, SceneEditor, MlabSceneModel

import Tkinter as tk
import tkFileDialog as fDog

class CrystalInteract(HasTraits):
    maxn  = 4
    #run_n = Range(low=0,high='maxn',value=0.1)
    #size  = Range(1,15,4)

    reloadB = Button('Load File')
    
    scene = Instance(MlabSceneModel, ())
    plot = Instance(PipelineBase)

    def __init__(self,c):
        HasTraits.__init__(self)
        self.crystal=c
        self.maxn = c.getNumberRuns()-1
        self.add_trait('run_n',Range(0,self.maxn,0))
        self.add_trait('size',Range(1,15,4))
        self.add_trait('scale',Range(.1,1.0,.5))
        #self.run_n.set(high=self.maxn)
        
    @on_trait_change('run_n,size,scale,scene.activated')
    def update_plot(self):
        if self.run_n > self.maxn:
            self.run_n = self.maxn
        x,y,z,c,sV=self.crystal.mkCrystalPoints(self.run_n,self.size)
        self.scene.mlab.clf()
        self.plot = self.scene.mlab.quiver3d(x, y, z, sV, sV, sV, scalars=c, mode="sphere", scale_factor=self.scale,vmin=0,colormap='jet')
        self.plot.glyph.color_mode = "color_by_scalar"
        self.plot.glyph.glyph_source.glyph_source.center = [0,0,0]
        

    @on_trait_change('reloadB')
    def reloadCrystal(self):
        root = tk.Tk()
        root.withdraw()
        fname = fDog.askopenfilename()
        self.crystal.loadFile(fname)
        self.maxn = self.crystal.getNumberRuns()-1
        #self.add_trait('run_n',Range(0,self.maxn,0))
        
    view = View(Item('scene', editor=SceneEditor(scene_class=MayaviScene),
                     height=250, width=300, show_label=False),
                Group('_', 'run_n', 'size','scale'),
                #Group('_', 'reloadB',),
                resizable=True,
                )


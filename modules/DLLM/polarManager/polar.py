# -*-mode: python; py-indent-offset: 4; tab-width: 8; coding: iso-8859-1 -*-
#  DLLM (non-linear Differentiated Lifting Line Model, open source software)
# 
#  Copyright (C) 2013-2015 Airbus Group SAS
# 
#  This program is free software; you can redistribute it and/or
#  modify it under the terms of the GNU General Public License
#  as published by the Free Software Foundation; either version 2
#  of the License, or (at your option) any later version.
# 
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
# 
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
# 
#  https://github.com/matthieu-meaux/DLLM.git
#
# @author : Francois Gallard
#
# - Local imports -
from OpenDACE.database.Database import Database
from differentiatedAeroShape import DifferentiatedAeroShape
from AeroElastAdj.kriging.kriging import Kriging
from numpy import matrix

class Polar(DifferentiatedAeroShape):

    def __init__(self,Sref,Lref,database=None,aeroShape=None,aOaList=None):
        '''
        Constructor for a polar.
        @param Sref : reference surface
        @param Lref : reference length        
        @param database : the dataBase in which the discrete values are stored, or None
        @param aeroShape : the aeroShape from which the polar should be built, or None if the database is already fulfilled
        @param aOaList : the angles of attack list on which  polar should be built, or None if the database is already fulfilled
        '''
        DifferentiatedAeroShape.__init__(self,Sref,Lref)
        self.__Sref=Sref
        self.__Lref=Lref
        self.__aeroShape=aeroShape
        self.__aOaList=aOaList
        self.__ClInterpolator=Kriging()
        self.__CdInterpolator=Kriging()
        self.__CmInterpolator=Kriging()

        if database is not None:
            self.__database=database
            self.__buildFromDatabase(database)
        elif aeroShape is not None and aOaList is not None:
            self.__database=Database()
            self.__buildFromAeroShape(aOaList,aeroShape)
        else:
            raise Exception, 'A polar initialization must be done either from a non-Null Database or a non-Null AeroShape and aOaList.'
            
    def __buildFromDatabase(self,database):
        '''
        Initializes the polar from a database in which are stored Cl, Cd, Cm for given aOa
        @param database : the database in which the data are stored
        '''
        self.__ClInterpolator.learn_from_database(database,'AoA','Cl')
        self.__CdInterpolator.learn_from_database(database,'AoA','Cd')
        self.__CmInterpolator.learn_from_database(database,'AoA','Cm')
                
    def __buildFromAeroShape(self,aOaList,aeroShape):
        '''
        Builds and stores a (aOa,Cl,Cd,Cm) polar based on the angles of attack list ans stores in the database
        @param aOaList : the list of angles of attack to build the polar on
        @param aeroShape : the aeroShape on which the polar is built
        '''
        for a in aOaList:
            self.__database.add_info(listin=['AoA'],AoA=a,Cd=aeroShape.Cd(a),Cl=aeroShape.Cl(a),Cm=aeroShape.Cm(a))
        
        self.__buildFromDatabase(self.__database)
        
    def Cl(self,alpha,beta=0.0,Mach=0.0):
        '''
        Lift coefficient function
        @param alpha: angle of Attack 
        @param type alpha : Float
        @param Mach : Mach number
        @param type Mach : Float
        @param beta : sideslip angle
        @param type beta : Float
        '''
        return self.__ClInterpolator.f(matrix([alpha]))[0,0]

    def Cd(self,alpha,beta=0.0,Mach=0.0):
        '''
        Drag coefficient function
        @param alpha: angle of Attack 
        @param type alpha : Float
        @param Mach : Mach number
        @param type Mach : Float
        @param beta : sideslip angle
        @param type beta : Float
        '''
        return (self.__CdInterpolator.f(matrix([alpha])))[0,0] 
    
    def Cm(self,alpha,beta=0.0,Mach=0.0):
        '''
        Pitch moment coefficient
        @param alpha: angle of Attack 
        @param type alpha : Float
        @param Mach : Mach number
        @param type Mach : Float
        @param beta : sideslip angle
        @param type beta : Float
        '''
        return (self.__CmInterpolator.f(matrix([alpha])))[0,0]
    
    def ClAlpha(self,alpha,beta=0.0,Mach=0.0):
        '''
        Sensibility of Cl to alpha.
        Done by finites differences by default.
        @param alpha: angle of Attack 
        @param type alpha : Float
        @param Mach : Mach number
        @param type Mach : Float
        @param beta : sideslip angle
        @param type beta : Float
        '''
        return self.__ClInterpolator.df(matrix([alpha]))[0,0]

    def CdAlpha(self,alpha,beta=0.0):
        '''
        Sensibility of Cd to alpha.
        Done by finites differences by default.
        @param alpha: angle of Attack 
        @param type alpha : Float
        @param Mach : Mach number
        @param type Mach : Float
        @param beta : sideslip angle
        @param type beta : Float
        '''
        return (self.__CdInterpolator.df(matrix([alpha])))[0,0] 
    
    def CmAlpha(self,alpha,beta=0.0):
        '''
        Sensibility of Cm to alpha.
        Done by finites differences by default.
        @param alpha: angle of Attack 
        @param type alpha : Float
        @param Mach : Mach number
        @param type Mach : Float
        @param beta : sideslip angle
        @param type beta : Float
        '''
        return (self.__CmInterpolator.df(matrix([alpha])))[0,0]
    
    def getDatabase(self):
        '''
        Accessor for the database in which the values are stored
        @return : the database
        '''
        return self.__database
    
    def export_to_XML(self,xmlFilePath):
        '''
        Exports the polar database to a XML file
        '''
        self.__database.export_to_XML(xmlFilePath)
        
    def get_aero_shape(self):
        """
        Accessor for the initial aeroShape.
        """
        return self.__aeroShape
    
    def get_AoA_list(self):
        """
        Accessor for the initial AoA list.
        """
        return self.__aOaList
        
    def exportPolarToFile(self,filePath,aOaList):
        '''
        Writes a polar AoA, Cl, Cd, LoD, Cm in a file
        @param file : the destination file
        @param aOaList : the angles of attack on which the polar is computes
        '''
        
        file=open(filePath,'w')
        for aOa in aOaList:
            Cl=self.Cl(aOa)
            Cd=self.Cd(aOa)
            Cm=self.Cm(aOa)
            LoD=Cl/Cd
            
            line = str(aOa)+" "+str(Cl)+" "+str(Cd)+" "+str(LoD)+" "+str(Cm)+"\n"
            file.write(line)
        file.close()
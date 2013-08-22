# -*-mode: python; py-indent-offset: 4; tab-width: 8; coding: iso-8859-1 -*-
# Copyright: Airbus

# - OpenDACE imports -
from Information import Information

# - XML imports -
from worms.XMLMod.XMLMod import PrettyPrint
from xml.dom.minidom import parse, getDOMImplementation
from xml.dom import Node

# - system import -
import os

class Database:
    """
    Class Database stores links between couples of informations.
    """

    BACKUP_FILE = 'backup_database.xml'
    
    def __init__(self, backup_file = None):
        """
        Class constructor.
        """
        self.__links        = []
        if backup_file is None:
            self.__backup_file = self.BACKUP_FILE
        else:
            self.__backup_file = backup_file

        # Try to load a previous database if regular save is asked
        # (Consider that if regular save is not asked, there's no need to restart from a crash...)
        if os.path.isfile(self.get_backup_file()):
            self.import_from_XML(self.get_backup_file())
    
    def get_backup_file(self):
        return self.__backup_file

    def __clean(self):
        self.__links = []
        
    def clean(self):
        self.__clean()

    def add(self, info1, info2, iprint=True):
        """
        Add a couple of informations.
        """
        # check items are type of Information
        if not isinstance(info1, Information):
            raise Exception, 'Database - add(): Error, wrong input type.'
        if not isinstance(info2, Information):
            raise Exception, 'Database - add(): Error, wrong output type.'

        # Check if couple already stored
        try:
            if (info1, info2) not in self.__links:
                # Add couple to the list of links
                self.__links.append((info1, info2))
        except:
            # Failed to compare (info may be complex objects)
            self.__links.append((info1, info2))
            if iprint:
                print 'Database - add: Warning, may duplicate objects stored.'

        return len(self.__links)

    def get_from_input(self, item_value, tag=Information.DEFAULT_TAG, eq_fun = None):
        """
        Returns a couple of informations from an input value.
        """
        # Go through all the links
        index=0
        for couple in self.__links:
            # Look for the desired value
            if tag in couple[0]:
                if (eq_fun is None and couple[0][tag] == item_value) or \
                   (eq_fun is not None and eq_fun(couple[0][tag], item_value)):
                    # Return couple and index (to index usage)
                    return couple, index
            index+=1

        return None, None
    
    def is_in_database(self, item_value, tag=Information.DEFAULT_TAG, eq_fun = None):
        c, i = self.get_from_input(item_value, tag=tag, eq_fun=eq_fun)
        return (c is not None)
            
    def add_info(self, iprint=True, listin=['x'], **kwargs):
        """
        Method to add info from variable arguments.
        """
        # Store input/output dictonnaries
        DictIn={}
        DictOut={}
        # Go through each argument
        for key, value in kwargs.items():
            # Check if arg is in or out (via listin)
            if key in listin:
                DictIn[key]=value
            else:
                DictOut[key]=value

        # Create Informations objects
        InfoIn=Information(DictIn)
        InfoOut=Information(DictOut)

        # Add link
        self.add(InfoIn,InfoOut, iprint)

    
    def get_from_index(self, index):
        """
        Returns a couple of informations from its index in the list.
        """
        return self.__links[index]
    

    def get_all_links(self):
        """
        Returns all the links stored in the database.
        """
        return self.__links
    
    def get_nb_links(self):
        """
        Returns the number of links
        """
        return len(self.__links)

    def get_all_in_outs(self,tagIn=None,tagOut=None):
        '''
        Returns a copy of the links stored in the database
        @param tagIn : the tag of input values to get
        @param tagOut : the tag of output values to get
        @return inputs : a list of the inputs
        @return outputs : a list of the outputs
        '''
        inputs=[]
        outputs=[]
        if(tagIn is None):
            for link in self.__links:
                inputs.append(link[0])
                outputs.append(link[1])
                
        if(tagIn is not None and tagOut is not None):
            for link in self.__links:
                inputs.append(link[0][tagIn])
                outputs.append(link[1][tagOut])
                
        if(tagIn is None and tagOut is not None):
            for link in self.__links:
                inputs.append(link[0])
                outputs.append(link[1][tagOut])
                
        if(tagIn is not None and tagOut is None):
            for link in self.__links:
                inputs.append(link[0][tagIn])
                outputs.append(link[1])             
            
        return inputs,outputs
    
    def get_all_x_y(self):
        '''
        Returns a copy of the values x & y stored in the database
        @return inputs : a list of the inputs matching x tag
        @return outputs : a list of the outputs values matching y tag
        '''        
        inputs=[]
        outputs=[]
        for link in self.__links:
            inputs.append(link[0]['x'])
            outputs.append(link[1]['y'].values())
        return inputs, outputs
            
    def export_to_XML(self, file_name):
        """
        Export the database to an XML file.
        """
        # Create xml document
        dom = getDOMImplementation()
        doc = dom.createDocument(None, 'Database', None)
        RootNode = doc.documentElement
        # Go through each link
        for ln in self.__links:
            # Create "link" node
            linkNode = doc.createElement('link')
            RootNode.appendChild(linkNode)
            # Export both Info part of the link
            ln[0].export_to_XML(doc, linkNode)
            ln[1].export_to_XML(doc, linkNode)

        # Write the xml file
        fid = open(file_name, 'w')
        PrettyPrint(doc,fid)
        fid.close()

    def import_from_XML(self, file_name, overwrite = True):
        """
        Import the database from an XML file.
        """
        # Clean the database if needed
        if overwrite:
            self.__clean()

        # Parse xml file
        dom = parse(file_name)
        RootNode = dom.childNodes[0]
        RootNodeName = RootNode.tagName
        if RootNodeName != 'Database':
            raise Excepetion, 'Database - import_from_XML: Error, given file is not a Database XML file.'

        # Go through each "link" node
        for LinkNode in RootNode.childNodes:
            if LinkNode.nodeType != Node.ELEMENT_NODE:
                pass
            else:
                # Store input/output informations
                info1 = None
                info2 = None
                for InfoNode in LinkNode.childNodes:
                    if InfoNode.nodeType != Node.ELEMENT_NODE:
                        pass
                    else:
                        # The two following lines may be replaced by : info = Information({})
                        # This code allows to store/load object derivated from Information
                        infoClass = eval(InfoNode.tagName)
                        info = infoClass({})
                        # Import the information
                        info.import_from_XML(InfoNode)
                        # Informations are ordered in the xml file
                        # So the first part of link node is input and the second output
                        if info1 is None:
                            info1 = info
                        else:
                            info2 = info

                # Check xml integrity and add link in database
                if info1 is None or info2 is None:
                    raise Exception, 'Database - import_from_XML: Error, corrupted XML file.'
                self.add(info1, info2, iprint=False)
    
    def var_equals(self, x1, x2):
        return (x1 == x2).all()
        






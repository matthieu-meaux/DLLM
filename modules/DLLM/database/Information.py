# -*-mode: python; py-indent-offset: 4; tab-width: 8; coding: iso-8859-1 -*-
# Copyright: Airbus

# - XML imports -
from xml.dom import Node

# - numpy import -
#import numpy

# - pickle import -
#import pickle


class Information(dict):
    """
    Class Information codes informations stored in the database.
    """

    DEFAULT_TAG = 'info'
            
    def __init__(self, data=None, tag=DEFAULT_TAG, **kwargs):
        """
        Class constructor.
        """
        # Check arguments
        if kwargs.keys()==[]:
            # Manage dictonnary argument
            if isinstance(data, dict):
                dict.__init__(self,data)
            else:
                # Manage single arg given
                dict.__init__(self)
                self[tag] = data
        else:
            # Manage multiple args given
            dict.__init__(self, kwargs)


    def export_to_XML(self, doc, parent_node):
        """
        Export the information to an XML file.
        """
        # Create the Information node
        InfoNode = doc.createElement('Information')
        parent_node.appendChild(InfoNode)
        # Go through each tag
        for key in self:
            # create the tag node
            keyNode = doc.createElement(key)
            InfoNode.appendChild(keyNode)
            # TODO: use pickle only for numpy?
            # Use pickle for every type allow to store only the value and no need for the type
            #val_str = pickle.dumps(self[key])
            val_str=str(self[key])
            TextNode = doc.createTextNode(val_str)
            keyNode.appendChild(TextNode)

        
    def import_from_XML(self, InfoNode):
        """
        Import the information from an XML file.
        """
        # Go through each child of the Information node (so each tag)
        for valNode in InfoNode.childNodes:
            if valNode.nodeType != Node.ELEMENT_NODE:
                pass
            else:
                # Get the value of the tag
                key = valNode.tagName
                val_str = str(valNode.childNodes[0].nodeValue)
                # TODO: use pickle for numpy only?
                # Will require to manage "None", string type and basic type
                #value = pickle.loads(val_str)
                try:
                    value=eval(val_str)
                except:
                    value=val_str
                self[key] = value

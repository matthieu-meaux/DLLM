import string

list_file=['mesh_ONERAM6_inv_Mach_03.dat','mesh_ONERAM6_inv_Mach_06.dat','mesh_ONERAM6_inv_Mach_08.dat']
list_corr=[0.0010251106,0.0008860518,0.0012233352]
for i,file_name in enumerate(list_file):
    fid=open(file_name,'r')
    lines=fid.readlines()
    fid.close()
    
    fid=open('corr_'+file_name,'w')
    for line in lines:
        words=string.split(line)
        if len(words)> 0:
            if words[0] != '#':
                words[2]=str(eval(words[2])-list_corr[i])
        new_line = string.join(words,' ')
        fid.write(new_line+'\n')
    fid.close()
        
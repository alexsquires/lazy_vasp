def cplaper( competition, elemental_references):
    
    ele_ref  = import_calculations_from_file(elemental_references)
    comp_ref = import_calculations_from_file(competition)
    #interest = import_calculations_from_file(material)
    
    elemental = cplaper_elements( elemental_references )
    
    competing = cplaper_competing( comp_ref, ele_ref )
    
    cplap_writer(energies, comp_ref, ele_ref)

    return elemental, competing
    
def cplap_writer(energies, competition, elements):
    
    cplap = []
    rad = []
    competing = []
    final = []
    res =[]

    
    for i in competing:
        j = (len (competition['{}'.format(i)].stoichiometry))
        final.append(j)
    

    for i in cplap:
        j = {v:k for k,v in i.items()}
        res.append(j)
    
    for key in ref_phas:
        m = str(key)
        competing.append(m)
        
    for i in competing:
        x = Calculation.scale_stoichiometry(competition['{}'.format(i)], 1)
        cplap.append(x)
    

    for i, j, k in zip( final, res , energies):
        rad.append(i)
        rad.append([j,k])

    with open('interim.dat', 'w') as file:
        for item in rad:
            file.write("%s\n" % item)
        

    f = open('interim.dat','r')
    filedata = f.read()
    f.close()

    newdata = filedata.replace("'","").replace("[","").replace("[","").replace("]","").replace(",","").replace("{","").replace("}","").replace(":","")


    f = open('cplap_input.dat','w')
    f.write(newdata)
    f.close()
    
def cplaper_elements( elemental_references ):
    
    ele_ref  = import_calculations_from_file(elemental_references)

    elements = []
    
    for key in ele_ref:                                                 
        n = str(key)
        elements.append(n)                                              

    elemental = []       
    

    for i in elements:                                                  
        n = float (ele_ref['{}'.format(i)].energy / ele_ref['{}'.format(i)].stoichiometry[i] )    
        y = {i:n}
        elemental.append(y)
    
    normalised_elemental_references = { k: v for d in elemental for k, v in d.items() }
    
    return normalised_elemental_references
    
def cplaper_competing( competing_phases, elemental_references ):
 
    energies = []
    myresults = []
    competing = []

    for key in ref_phas:
        m = str(key)
        competing.append(m) 
    
    
    for compound in competing:
        for element in elemental:
            h =  ref_phas['{}'.format(compound)].stoichiometry['{}'.format(element)] * elemental['{}'.format(element)] 
            myresults.append(h)
            it = iter(myresults)
            k = [a + b + c + d for a,b,c,d in zip(it,it,it,it)]
            
    for compound,i in zip (competing,k):
         
        result = ( ( ref_phas['{}'.format(compound)].energy) - i ) / sum(ref_phas['{}'.format(compound)].stoichiometry.values())
            
        energies.append(result)    
    
          
    
    return energies
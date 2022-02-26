def menu():
    '''
    The function displays a menu for the available options for the operations available for the PDB file.
    '''
    menu_dict = {1 : '*' * 80, 
            2 : ' PDB FILE ANALYZER',
            3 : '*' * 80,
            4 : ' Select an option from below:',
            5 : '',
            6 : '     1) Open a PDB File                      (O)',
            7 : '     2) Information                          (I)',
            8 : '     3) Show histogram of amino acids        (H)',
            9 : '     4) Display Secondary Structure          (S)',
            10: '     5) Exit                                 (Q)',
            11: '',
            12: 'Current PDB: None '.rjust(80),
            13: '*' * 80}
    
    for line in menu_dict.values():
            print('*%s*' % line.ljust(80))
    
    
#Verification of file before loading functions 
def pdb_verification1(pdbfile):
    '''
    Verifies the file follows the pdb format using the length of columns and first 6 columns of the first line.
    '''
    first_line = pdbfile.readline()
    first_line = first_line.replace('\n' , '')
    verify = len(first_line) == 80 and first_line.startswith('HEADER')
    return verify

def pdb_verification2(pdbfile):
    '''
    Verifies the length of each line and of the maximum limit of the first record and returns a boolean.
    '''
    pdbfile = pdbfile.readlines()
    for line in pdbfile:
        record_name = line.split()[0]
        verify2=len(line)==81 and len(record_name) <= 6
    return verify2

#Function for verifying and opening and loading the file to memory
def open_file(path):
    '''
    Loads the file into memory after veryfying it is a pdb file.
    '''

    try:
        with open(path , 'r') as pdbfile:
            if path.endswith('.pdb') and pdb_verification1(pdbfile) and pdb_verification2(pdbfile):
                validpath = path
                return validpath
            else:
                print('Please provide a valid pdbfile')
                
    except FileNotFoundError:
        print('Path does not exist')
    
    
#Function for extracting the file name to be inserted in the modified menu

def filename(r):
    '''
    Extracts filename to be displayed on the menu after loading file into memory.
    '''
    try:
        with open(r, 'r') as pdbfile2:
            line1 = pdbfile2.readline()
            name = line1[62:66] + '.pdb'
            return name
    except:
        pass
    

#Modified menu displaying the loaded file.The menu should be displayed for selection of any options after a file has been loaded.It is displayed after execution of any option selected. 
        
def menu_mod(name):
    '''
    The function displays a menu for the operations available on the PDB file with a display of the loaded filename.
    '''
    
    menu_dict = {1 : '*' * 80, 
            2 : ' PDB FILE ANALYZER',
            3 : '*' * 80,
            4 : ' Select an option from below:',
            5 : '',
            6 : '     1) Open a PDB File                      (O)',
            7 : '     2) Information                          (I)',
            8 : '     3) Show histogram of amino acids        (H)',
            9 : '     4) Display Secondary Structure          (S)',
            10: '     5) Exit                                 (Q)',
            11: '',
            12: 'Current PDB: None '.replace('None', name).rjust(80),
            13: '*' * 80}
    
    for line in menu_dict.values():
        print('*%s*' % line.ljust(80))
                                            
    
#OPTION2
#Obtains the formatted string of chains in the sequence
def chain_list(chain_id_list):
    '''
    The function carries out formatting of itens from one list returns a formatted list 
    '''
    pr_chain_list = '' 
    if len(chain_id_list) == 1:
        pr_chain_list+=chain_id_list[0]
    elif len(chain_id_list) == 2:
        pr_chain_list+=chain_id_list[0] + ' ' + 'and' + ' ' + chain_id_list[1]
    else:
        for i in chain_id_list:
            if chain_id_list.index(i) != len(chain_id_list)-1 and  chain_id_list.index(i) != len(chain_id_list)-2:
                pr_chain_list+= i+ ',' + ' ' 
        pr_chain_list+= chain_id_list[-2] + ' ' + 'and' + ' ' + chain_id_list[-1]
    return pr_chain_list


#Function to extract and output information in a specified format.
def option_2():
    '''
    The function follows a series of steps that extracts various columns from a file and stores them in lists used to make accessible the information.
    '''
    x = 0
    while x< len(chain_id_list):
        numbers = []
        for line in infofile_list:
            if line.startswith('SEQRES'):
                if line[11] == chain_id_list[x]:
                    number = line[14:17]
                    if number in numbers:
                        continue

                    else:
                        numbers.append(number)



        helix= []
        try:
            for line in infofile_list:
                if line.startswith('HELIX'):
                    if line[19] == chain_id_list[x]:
                        helix_numbers = line[9]
                        helix.append(helix_numbers)
            helix = len(helix)
        except ValueError:
            helix = 0



        count = 0
        sheets = []
        for line in infofile_list:
            if line.startswith('SHEET'):
                if line[21] == chain_id_list[x]:
                    if line[17:21] in sheets:
                        continue
                    else:
                        sheets.append(line[17:20])
                        count = count+1


        codons_dict = {
            'ALA' : 'A' , 'ARG' : 'R' , 'ASN' : 'N' , 
            'ASP' : 'D' , 'CYS' : 'C' , 'GLU' : 'E' ,
            'GLN' : 'Q' , 'GLY' : 'G' , 'HIS' : 'H' ,
            'ILE' : 'I' , 'LEU' : 'L' , 'LYS' : 'K' ,
            'MET' : 'M' , 'PHE' : 'F' , 'PRO' : 'P' ,
            'SER' : 'S' , 'THR' : 'T' , 'TRP' : 'W' ,
            'TYR' : 'Y' , 'VAL' : 'V'}


        codons_string = ''             
        for line in infofile_list:
            if line.startswith('SEQRES'):
                if line[11] == chain_id_list[x]:
                    codons_string+=line[19:]
       

        codons_string =codons_string.replace('  ', '').replace('\n', ' ')
        codons_list = codons_string.split()


        sequence_list = []
        for residue in codons_list:
            if len(residue) == 3 and residue in codons_dict:
                sequence_list.append(codons_dict[residue])
            elif len(residue) == 1:
                sequence_list.append(residue)
        sequence_str = ''.join(sequence_list)
        
        sequence = textwrap.fill(sequence_str,width=50) 
        sequence_wr = textwrap.indent(sequence,predicate=lambda sequence_str: not sequence.splitlines()[0] in sequence_str,prefix='\t       ') 
        
        
        print(f' - Chain {chain_id_list[x]}')
        
        print('%4sNumber of aminoacids:%7s' %(' ',numbers[0]))
        print('%4sNumber of helix:%12s' %(' ',helix))
        print('%4sNumber of sheet:%12s' %(' ',count))
        print(f'    Sequence:  {sequence_wr}')
        

        x = x+1

#OPTION3
#creates sorted dictionaries used to generate histograms 
def option3(aa_count_dict):
        '''
        Sorts the dictionary containing amino acids in the described order and represents the number as a histogram 
        '''
        if order_choice == 'an':
            numeric_asc_dict = dict(sorted(aa_count_dict.items(),key =lambda item:item[1]))
            for item in numeric_asc_dict:
                hist_asc_count = '*' * numeric_asc_dict[item]
                print('%s (%3d) : %s' % (item, numeric_asc_dict[item], hist_asc_count))

        elif order_choice == 'dn':
            numeric_desc_dict = dict(sorted(aa_count_dict.items(),key = lambda item:item[1],reverse= True))
            for item in numeric_desc_dict:
                hist_desc_count = '*' * numeric_desc_dict[item]
                print('%s (%3d) : %s' % (item, numeric_desc_dict[item], hist_desc_count))


        elif order_choice == 'aa':
            alphabetical_asc_keys = sorted(aa_count_dict.keys())
            alphabetical_asc_dict = {key:aa_count_dict[key] for key in alphabetical_asc_keys}
            for item in alphabetical_asc_dict:
                hist_aasc_count = '*' * alphabetical_asc_dict[item] 
                print('%s (%3d) : %s' % (item, alphabetical_asc_dict[item], hist_aasc_count))


        elif order_choice == 'da':
            alphabetical_desc_keys = sorted(aa_count_dict.keys(), reverse=True)
            alphabetical_desc_dict = {key:aa_count_dict[key] for key in alphabetical_desc_keys}
            for item in alphabetical_desc_dict:
                hist_desc_count = '*' * alphabetical_desc_dict[item] 
                print('%s (%3d) : %s' % (item, alphabetical_desc_dict[item], hist_desc_count))
                
#OPTION4
def option_4(chain_id_list, infofile_list):
    '''
    Utilizes a translated sequence and special characters to represent the secondary structure of a protein chain.
    '''
    x = 0 
    while x< len(chain_id_list):
        
        
        codons_string = '' 
        
        for line in infofile_list:
            if line.startswith('SEQRES'):
                if line[11] == chain_id_list[x]:
                    codons_string+=line[19:]


        codons_string =codons_string.replace('  ', '').replace('\n', ' ')
        codons_list = codons_string.split()
        
        codons_dict = {
            'ALA' : 'A' , 'ARG' : 'R' , 'ASN' : 'N' , 
            'ASP' : 'D' , 'CYS' : 'C' , 'GLU' : 'E' ,
            'GLN' : 'Q' , 'GLY' : 'G' , 'HIS' : 'H' ,
            'ILE' : 'I' , 'LEU' : 'L' , 'LYS' : 'K' ,
            'MET' : 'M' , 'PHE' : 'F' , 'PRO' : 'P' ,
            'SER' : 'S' , 'THR' : 'T' , 'TRP' : 'W' ,
            'TYR' : 'Y' , 'VAL' : 'V'
        }

        
        sequence_list = []
        for residue in codons_list:
            if len(residue) == 3 and residue in codons_dict:
                sequence_list.append(codons_dict[residue])
            elif len(residue) == 1:
                sequence_list.append(residue)
        sequence_str = ''.join(sequence_list)


        helix_indices= []
        for line in infofile_list:
            if line.startswith('HELIX'):
                if line[19] == chain_id_list[x]:
                    helix_index1=line[22:25]
                    helix_index2= line[34:37]
                    helix_indices.append([(int(helix_index1)-1),int(helix_index2)])
        

        helix_secondary_string = ''
        helix_secondary_string+=sequence_str
        hlist=helix_secondary_string
        hlist = list(hlist)
                     

        for index in helix_indices:
            for i in range(index[0],index[1]):
                hlist[i] = '/'
                     

        sheet_indices = []
        for line in infofile_list:
                if line.startswith('SHEET'):
                    if line[32] == chain_id_list[x]:
                        sheet_index1 = line[23:26]
                        sheet_index2 = line[34:37]
                        sheet_indices.append([(int(sheet_index1)-1),int(sheet_index2)])
                     
        
        for index in sheet_indices:
            for i in range(index[0],index[1]):
                hlist[i] = '|'
        
        
        secondary_structure = ''.join(hlist)
        for char in secondary_structure:
            if char in sequence_str:
                secondary_structure = secondary_structure.replace(char, '-').replace(',', '')

   
        sheet_strandnumber_list = []
        for line in infofile_list:
                if line.startswith('SHEET'):
                    if line[21] == chain_id_list[x]:
                        sheet_strandnumber_list.append(line[7:10])
        

        sheet_identifier_list = []
        for line in infofile_list:
                if line.startswith('SHEET'):
                    if line[21] == chain_id_list[x]:
                        sheet_identifier_list.append(line[11:14])
        print(sheet_identifier_list)

        sheet_index1_list = []
        for line in infofile_list:
                if line.startswith('SHEET'):
                    if line[21] == chain_id_list[x]:
                        sheet_index1 = line[23:26]
                        sheet_index1_list.append(int(sheet_index1)-1)       

        helix_identifier_list= []
        for line in infofile_list:
            if line.startswith('HELIX'):
                if line[19] == chain_id_list[x]:
                    helix_identifier_list.append(line[11:14])
        print(helix_identifier_list)


        helix_index1_list = []
        for line in infofile_list:
            if line.startswith('HELIX'):
                if line[19] == chain_id_list[x]:
                    helix_index1 = line[22:25]
                    helix_index1_list.append(int(helix_index1)-1)


        tag_list = list(secondary_structure)[:]
        
        
        for i,j,k in zip(sheet_strandnumber_list,sheet_identifier_list, sheet_index1_list):
                print(i,j)
                if len(i.strip()) == 1 and len(j.strip()) == 1:
                    tag_list[int(k)] = i.strip()
                    tag_list[int(k)+1] = j.strip()
                    

        index_list = []
        for i,j in zip(helix_identifier_list,helix_index1_list):
            if len(i.strip()) == 1:
                tag_list[int(j)] = i.strip()
            elif len(i.strip()) == 2:
                index_list.append(list(i.strip()))
                tag_list[int(j)] = index_list[0][0]
                tag_list[int(j)+1] = index_list[0][1]
            elif len(i.strip()) == 3:
                index_list.append(list(i.strip()))
                tag_list[int(j)] = index_list[0][0]
                tag_list[int(j)+1] = index_list[0][1]           
            
            print(tag_list)
                            
        tag_string = ''.join(tag_list)
        for char in tag_string:
            if char == '/' or char == '|' or char == '-':
                tag_string = tag_string.replace(char, ' ').replace(',', '')

           
        print(f'Chain {chain_id_list[x]}:')
        print('(1)')
        for i in range(0,len(sequence_str),80):
            print(sequence_str[i:i+80]+'\n'+secondary_structure[i:i+80]+'\n'+tag_string[i:i+80]+'\n')
        print('\n')
        print(f'({len(secondary_structure)})')
        print('\n')
        x=x+1
        



#""""""""""""""""""""""""""""""""""""""""""""""""""""""""""PROGRAMCODE"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
import textwrap 

#""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
#"""""""""""""""""""""""""""""""""""""""""""""""""""""""OPTION1(main_menu)""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
#The menu is displayed for the user to select an option to load the file. The user is restricted to choosing option 1 or 'O' or q to exit the program.The user can only access other options after loading a file. A valid path to the file MUST be provided and is checked by a function that verifies the path provided by user exists and proceeds to check for a '.pdb' extension. The function also verifies that the file loaded is not empty and follows the pdb format.
#"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
 

while True:
    
    menu()
    print('Select 1 or O to load a file or Q to exit')
    menu_choice = input(': ')
    options = ['1', 'O'.upper(),'Q'.upper()]
    
    if menu_choice.upper() == 'Q':
        break
    
    if menu_choice.upper() not in options: 
        print('Invalid choice.')
        continue
            

    if menu_choice.upper() == '1' or menu_choice.upper() == 'O':
        print(f': {menu_choice.upper()}')
        path = input('Enter a valid path for a PDB file: ')
        if open_file(path) != path:
                print('Invalid choice')
                continue
            
                
        else:
            c = open_file(path)
            name = filename(c)
            print(f'The file {name} has been successfuly loaded ')

            

#""""""""""""""""""""""""""""""""""""""""""""OPTION1(modified_menu)"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
# A modified menu displaying the name of the loaded file is displayed. The user can load another file after confirming to replace the already loaded file by reselecting the option 1 or   'O'. 
#"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

    while True:
        try:
            name = filename(c)
            menu_mod(name) 

        except:
            pass
        
        menu_mod_choice = input(': ')

        valid_options = ['1','2','3','4','5','O','I','H','S','Q']
        if menu_mod_choice.upper() not in valid_options:
            print('Invalid choice')
            continue

        if menu_mod_choice =='1' or menu_mod_choice.upper() =='O':
            rep_options= ['yes', 'no', 'q']
            rep_choice = input('''Are you sure you want to replace the current file?
                           Select :
                               yes -  to replace the current file 
                               no  -  to retain the current file
                               q   -  to return to menu
                            : ''')

            if rep_choice not in rep_options:
                print('Invalid choice')
                continue            

            if rep_choice.lower() == 'yes':
                path = input('Enter a valid path for a PDB file: ')

                while open_file(path) != path:
                    print('Input valid path or q to return to main menu')
                    path = input('Enter a valid path for a PDB file: ')  
                    if path == 'q'.lower():
                            break
                else:
                    c = open_file(path) 
                    name = filename(c)



            elif rep_choice.lower() == 'no':
                print('Proceed with current file. Select option from menu below')
            else:
                continue


#""""""""""""""""""""""""""""""""""""""""""""""""OPTION2(Information)"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
# The user is provided with information on the file loaded. The program considers each column as a record and extracts them through indexing and presents the information to the user.
#Information on title(len 80 per line),the chains in the file,number of aminoacids(width7),number of helices and sheets(width 12) in a chain and the sequence of the chain(len 50 per line) in one letter amino acid code is presented.  
#"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

        elif menu_mod_choice == '2' or menu_mod_choice.upper() == 'I':
            print(f': {menu_mod_choice.upper()}')
            print(f'PDB File: {filename(c)}')

            with open(path) as infofile:
                infotitle = ''
                infofile_list = infofile.readlines()
                for line in infofile_list:
                    if line.startswith('TITLE'):
                        infotitle+=line
                infofile_title =infotitle.replace('  ', '').replace('TITLE', '').replace('\n', '')

                infofile_title_wrap = textwrap.fill(infofile_title,width=80) 
                print(f'Title: {infofile_title_wrap}')


            chain_id_list = ''
            for line in infofile_list:
                    if line.startswith('SEQRES'):
                        chain = line[11]
                        if chain in chain_id_list:
                            continue
                        else:
                            chain_id_list = chain_id_list+chain


            print(f'CHAINS: {chain_list(chain_id_list)}')


            option_2()

#""""""""""""""""""""""""""""""""""""""""""""""""OPTION3(Histogram)""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
# This option allows the user to view the frequency of an amino acid in the sequence,irrespective of the chains. The user can choose to order the histogram by different methods. The #number is justified within(width 3).
#""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

        elif menu_mod_choice == '3' or  menu_mod_choice.upper() == 'H':
            quit=''

            while quit !='q':

                print(': %s' %(menu_mod_choice.upper()))      
                # order_choice = input('Choose an option to order by: ')           
                print('%3snumber of amino acids - ascending (an)' %' ')
                print('%3snumber of aminoacids - descending (dn)' %' ')
                print('%3salphabetically - ascending        (aa)' %' ')
                print('%3salphabetically - descending       (da)' %' ')
                order_choice = input('Choose an option to order by or q to quit the submenu: ') 
                print(f'order by: {order_choice.lower()}')
                print('\n')


                order_options = ['aa', 'an', 'dn', 'da','q']
                if order_choice  not in order_options:
                    print('Invalid choice')
                elif order_choice=='q':
                    quit+=order_choice

                else:

                    codons_string = '' 
                    with open(path) as infofile:
                        infofile_list = infofile.readlines()
                        for line in infofile_list:
                            if line.startswith('SEQRES'):
                                    codons_string+=line[19:]

                    codons_string =codons_string.replace('  ', '').replace('\n', ' ')
                    codons_list = codons_string.split()
                    codons_set =set(codons_list)

                    count_list = []
                    for amino_acid in codons_set:
                        count = codons_list.count(amino_acid)
                        count_list.append(count)

                    aa_count_dict= dict((aa,number) for aa,number in zip(codons_set,count_list))
                    option3(aa_count_dict)


        
#""""""""""""""""""""""""""""""""""""""""OPTION4(Secondarystructure)""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
#Provides the user with an indication of the secondary structure components and their locations. The structure is formed by helices and sheets of varying lengths. Corresponding identifier tags are also aligned below their respective substructures.
#"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

        elif menu_mod_choice == '4' or menu_mod_choice.upper() == 'S':
            print(': %s' %menu_mod_choice.upper())
            print('Secondary structure of the PDB id %s:' %filename(c))

            chain_id_list = ''
            with open(path) as infofile:
                infofile_list = infofile.readlines()
            for line in infofile_list:
                    if line.startswith('SEQRES'):
                        chain = line[11]
                        if chain in chain_id_list:
                            continue
                        else:
                            chain_id_list = chain_id_list+chain

            option_4(chain_id_list,infofile_list)

#"""""""""""""""""""""""""""""""""""""""""""""""""""""""OPTION5(exit)"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
# Allows the user to exit the program or return to the menu
#"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
        
        elif menu_mod_choice == '5' or menu_mod_choice.upper() == 'Q':
            exit_choice = input('Do you want to exit(E) or do you want go back to the menu (M): ')
            print(f': {exit_choice.upper()}')
            if exit_choice.upper() == 'M':
                continue
            elif exit_choice.upper() == 'E':
                break

#"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""




            
    
        
        
        



    



        
        







    

    

                                             
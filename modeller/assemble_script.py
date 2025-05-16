    """
    Script to force assembly of permuted protein. Written by Jaka Snoj (jaka.snoj@ki.si)
    """

from modeller import *
from modeller.scripts import complete_pdb
from modeller.optimizers import conjugate_gradients, molecular_dynamics, quasi_newton, actions
from modeller.automodel import *
import random as rd
env = environ()
env.libs.topology.read('${LIB}/top_heav.lib') 
env.libs.parameters.read('${LIB}/par.lib')
env.edat.dynamic_sphere = True
scale = 1 
for i in range(0,1000):
    for name in ['SplitLuc_1', 'SplitLuc_2', 'LCB3_1', 'LCB3_2']:
        mdl=complete_pdb(env,'{}.pdb'.format(name),transfer_res_num=True)
        atmsel = selection(mdl)
        x = scale*rd.uniform(-1,1)
        y = scale*rd.uniform(-1,1)
        z = scale*rd.uniform(-1,1)
        print('x,y,z',x,y,z)
        atmsel.translate([x,y,z] )
       
        mdl.write(file='{}_{}.pdb'.format(name,i))
        mdl.write(file='{}_current.pdb'.format(name))
 
    
#VtgyrlFeeil in MlfrvtIn (V in I, F in M)
    SplitLuc_1_S_N = 73 #V
    SplitLuc_1_L_N = 496 #I
    
    SplitLuc_1_S_C = 80 #F
    SplitLuc_1_L_C = 490 #M
    
    SplitLuc_2_S_N = 327 #V
    SplitLuc_2_L_N = 242 #I
    
    SplitLuc_2_S_C = 335 #F
    SplitLuc_2_L_C = 236 #M
    

    for model in range(1,2):
        class AlphaModel(automodel):
            def special_restraints(self, aln):
                rsr = self.restraints
                at = self.atoms
                                                               
                rsr.add(forms.gaussian(group=physical.xy_distance,     
                                   feature=features.distance(self.residues[SplitLuc_1_S_N].atoms['CA'],
                                                            self.residues[SplitLuc_1_L_N].atoms['CA']), 
                                   mean=float(5.5), stdev=float(1)))
                                   
                rsr.add(forms.gaussian(group=physical.xy_distance,     
                                   feature=features.distance(self.residues[SplitLuc_1_S_C].atoms['CA'],
                                                            self.residues[SplitLuc_1_L_C].atoms['CA']), 
                                   mean=float(5.5), stdev=float(1)))
                rsr.add(forms.gaussian(group=physical.xy_distance,     
                                   feature=features.distance(self.residues[SplitLuc_2_S_N].atoms['CA'],
                                                            self.residues[SplitLuc_2_L_N].atoms['CA']), 
                                   mean=float(5.5), stdev=float(1)))
                                   
                rsr.add(forms.gaussian(group=physical.xy_distance,     
                                   feature=features.distance(self.residues[SplitLuc_2_S_C].atoms['CA'],
                                                            self.residues[SplitLuc_2_L_C].atoms['CA']), 
                                   mean=float(5.5), stdev=float(1)))                                                                                                 
                                   
                rsr.add(forms.gaussian(group=physical.xy_distance,     
                                   feature=features.distance(self.residues[SplitLuc_2_S_C].atoms['CA'],
                                                            self.residues[SplitLuc_2_S_N].atoms['CA']), 
                                   mean=float(19.7), stdev=float(1)))
                rsr.add(forms.gaussian(group=physical.xy_distance,     
                                   feature=features.distance(self.residues[SplitLuc_1_S_C].atoms['CA'],
                                                            self.residues[SplitLuc_1_S_N].atoms['CA']), 
                                   mean=float(19.7), stdev=float(1))) 
                                                                                                
                rsr.add(secondary_structure.strand(self.residue_range(SplitLuc_1_S_N,SplitLuc_1_S_C)))
                rsr.add(secondary_structure.strand(self.residue_range(SplitLuc_2_S_N,SplitLuc_2_S_C)))
                        
            def get_model_filename(self, sequence, id1, id2, file_ext):
                print(self, sequence, id1, id2, file_ext)

                return "pLUC_1_slow_close_betaRestr_-{id2:02}-{id1}{file_ext}".format(id2=i,id1=model,file_ext=file_ext) 

        a = AlphaModel(env, alnfile='PLUC.ali', 
              knowns=('SplitLuc_1_S', 'SplitLuc_1_L', 'LCB3_1', 'SplitLuc_2_S', 'SplitLuc_2_L', 'LCB3_2'), sequence='pluc',
                      assess_methods=[assess.DOPE,
                                      #soap_protein_od.Scorer(),
                                      assess.GA341])
        a.md_level = refine.slow  
        a.starting_model = 1
        a.ending_model = 1
        a.max_molpdf 
        a.make()
    
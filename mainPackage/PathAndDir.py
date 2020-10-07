import os, sys
from pathlib import Path
#import pathlib
from GraphPackage.Graph_Config import CutOffContact
#Dir2Base = '../MyDataSet/DynDom/Perfect/BackUp/WithoutScaler_10.5'

Root_Dir = Path(__file__).parent.parent.absolute()
Dir2Base = os.path.join(str(Root_Dir),'Output') #'../Output'
#Dir2Base = os.path.join(Dir2Base, str(CutOffContact))

#Dir2SelectedDynDomEntry     = '../MyDataSet/DynDom/SelectedDynDomEntry'
#Path2SelectedSerial         = '../MyDataSet/DynDom/SerialList.txt'
#Dir2SelectedEntryFinal      = '../MyDataSet/DynDom/SelectedEntryFinal'

Path2ViterbiJar = os.path.join(str(Root_Dir),'Script/ViterbiJar/ViterbiAlgorithm.jar')
Dir2ClusterGraph            = os.path.join(Dir2Base, 'ClusterGraph') #'../MyDataSet/DynDom/Perfect/Graph'
Dir2LineGraph               = os.path.join(Dir2Base, 'LineGraph') #'../MyDataSet/DynDom/Perfect/LineGraph'
Dir2ViterbiFeature          = os.path.join(Dir2Base, 'ViterbiFeature') #'../MyDataSet/DynDom/Perfect/ViterbiFeature'
Dir2ViterbiOutFile          = os.path.join(Dir2Base, 'ViterbiOutFile') #'../MyDataSet/DynDom/Perfect/ViterbiOutFile'
Dir2GraphFigure             = os.path.join(Dir2Base, 'Figure') #'../MyDataSet/DynDom/Perfect/Figure'
Dir2Figure_Result           = os.path.join(Dir2Base, 'Figure_Result')
Dir2Figure_Result_Trivial   = os.path.join(Dir2Base, 'Figure_Result_Trivial')
Dir2Result                  = os.path.join(Dir2Base, 'Result') #'../MyDataSet/DynDom/Perfect/Result'
Dir2FeatureHistogram        = os.path.join(Dir2Base, 'Histogram')
Dir2G_Feature               = os.path.join(Dir2Base, 'G_Feature')
Dir2TmpFile                 = os.path.join(Dir2Base, 'TmpFile')
Dir2RMSD_Analysis           = os.path.join(Dir2Base, 'RMSD_Analysis')
Dir2Result_Cluster          = os.path.join(Dir2Base, 'Result_Cluster')
Dir2Membership              = os.path.join(Dir2Base, 'Membership')
Dir2ProteinGraph            = os.path.join(Dir2Base, 'ProteinGraph')


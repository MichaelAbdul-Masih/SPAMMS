
def line_bounds():
    x = {'HEPS':[3950, 3990], 'HEI4026':[4008, 4040], 'NIV4058':[4053, 4062], 'CIII4069':[4058, 4078], 'HDELTA':[4080, 4130], 'HEI4121':[4105, 4130], 'HEI4143':[4135, 4155], 'CIII4187':[4175, 4195], 'HEII4200':[4182, 4210],
    'HGAMMA':[4323, 4359], 'NIII4379':[4370, 4397], 'HEI4387':[4370, 440], 'DIB':[4415.0, 4454.4], 'HEI4471':[4467, 4475], 'NIIIqua':[4507, 4524], 'HEII4541':[4534, 4548], 'NIIItrip':[4626, 4645], 'CIII4650':[4638, 4661],
    'HEII4686':[4679, 4694], 'HEI4713':[4703, 4723], 'HBETA':[4847, 4877], 'HEI4922':[4915, 4930], 'HEI5016':[5008, 5024], 'HEII5411':[5400, 5420], 'OIII5592':[5584, 5602], 'CIV5801':[5790, 5824],
    'HEI5875':[5866, 5885], 'HALPHA':[6545, 6580], 'HEI6678':[6670, 6695], 'HEI7065':[7050, 7080], 'HEI170':[16975, 17040], 'HEI212':[21110, 21140], 'HEII169':[16900, 16960], 'HEII218':[21860, 21925]}
    return x

def abundance_dictionary():
    x = {
    'He':['HEI4026', 'HEI4143', 'HEI4387', 'HEI4471', 'HEI4713', 'HEI4922', 'HEI5016', 'HEI5048', 'HEI5875', 'HEI6678', 'HEI7065', 'HEI170', 'HEI212', 'HEII1640', 'HEII4200', 'HEII4541', 'HEII4686', 'HEII6406', 'HEII6527', 'HEII6683', 'HEII169', 'HEII218', 'HALPHA', 'HBETA', 'HGAMMA', 'HDELTA', 'HEPS'],
    'C':['CIII1176', 'CIII1620', 'CIII4069', 'CIII5695', 'CIII5697', 'CIV1169', 'CIV1512', 'CIV1548', 'CIV5801'],
    'N':['NII3995', 'NII4447', 'NII4601', 'NII4607', 'NII4621', 'NII4630', 'NII4643', 'NII5005', 'NII5007', 'NIII1750', 'NIII4003', 'NIII4379', 'NIII4511', 'NIII4515', 'NIII4518', 'NIII4523', 'NIII4527', 'NIII4531', 'NIII4535', 'NIII4547', 'NIII4634', 'NIII4640', 'NIII4641', 'NIIIqua', 'NIIItrip', 'NIV1718', 'NIV3480', 'NIV4058', 'NIV6380', 'NV1250', 'NV4603', 'NV4619', 'NVall', 'HDELTA', 'HEII4200'],
    'O':['OIII5592', 'OIV1340', 'OV1371', 'OVI1031', 'OVI1037']
    }
    return x

full_FASTWIND_line_list = ['CIII1176', 'CIII1620', 'CIII4069', 'CIII5695', 'CIII5697', 'CIV1169', 'CIV1512', 'CIV1548', 'CIV5801', 'HALPHA', 'HBETA', 'HDELTA', 'HEI4026', 'HEI4143', 'HEI4387', 'HEI4471', 'HEI4713', 'HEI4922', 'HEI5016', 'HEI5048', 'HEI5875', 'HEI6678', 'HEI7065', 'HEI170', 'HEI212', 'HEII1640', 'HEII4200', 'HEII4541', 'HEII4686', 'HEII6406', 'HEII6527', 'HEII6683', 'HEII169', 'HEII218', 'HEPS', 'HGAMMA', 'NII3995', 'NII4447', 'NII4601', 'NII4607', 'NII4621', 'NII4630', 'NII4643', 'NII5005', 'NII5007', 'NIII1750', 'NIII4003', 'NIII4379', 'NIII4511', 'NIII4515', 'NIII4518', 'NIII4523', 'NIII4527', 'NIII4531', 'NIII4535', 'NIII4547', 'NIII4634', 'NIII4640', 'NIII4641', 'NIIIqua', 'NIIItrip', 'NIV1718', 'NIV3480', 'NIV4058', 'NIV6380', 'NV1250', 'NV4603', 'NV4619', 'NVall', 'OIII5592', 'OIV1340', 'OV1371', 'OVI1031', 'OVI1037', 'PV1128a', 'PV1128b', 'SIIV1393', 'SIIV4654', 'SIIV6667', 'SIIV6701']

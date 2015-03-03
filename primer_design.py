import primer3
import sys

# example usage: python primer_design.py /Users/eric/d/j/cureffilab/media/2015/03/sgrna_designs_with_context.txt > ~/d/j/cureffilab/media/2015/03/sgrna_designs_with_primers.txt

f = open(sys.argv[1],mode='r')

primer3_general_settings =  {
        'PRIMER_OPT_SIZE': 20,
        'PRIMER_PICK_INTERNAL_OLIGO': 1,
        'PRIMER_INTERNAL_MAX_SELF_END': 8,
        'PRIMER_MIN_SIZE': 18,
        'PRIMER_MAX_SIZE': 36,
        'PRIMER_OPT_TM': 60.0,
        'PRIMER_MIN_TM': 57.0,
        'PRIMER_MAX_TM': 63.0,
        'PRIMER_MIN_GC': 20.0,
        'PRIMER_MAX_GC': 80.0,
        'PRIMER_MAX_POLY_X': 100,
        'PRIMER_INTERNAL_MAX_POLY_X': 100,
        'PRIMER_SALT_MONOVALENT': 50.0,
        'PRIMER_DNA_CONC': 50.0,
        'PRIMER_MAX_NS_ACCEPTED': 0,
        'PRIMER_MAX_SELF_ANY': 12,
        'PRIMER_MAX_SELF_END': 8,
        'PRIMER_PAIR_MAX_COMPL_ANY': 12,
        'PRIMER_PAIR_MAX_COMPL_END': 8,
        'PRIMER_PRODUCT_SIZE_RANGE': [[150,200]],
        'PRIMER_MAX_NS_ACCEPTED': 1
    }

header = f.readline().strip('\n').split('\t')
print '\t'.join(header) + '\t' + '\t'.join(['left','right','left_tm','right_tm','product_size'])
for line in f.readlines():
    data = dict(zip(header,line.strip('\n').split('\t')))
    sequence_template = data['seq_context']
    sequence_id = data['symbol'] + '_' + data['spacer']
    specs = {'SEQUENCE_ID': sequence_id, 'SEQUENCE_TEMPLATE': sequence_template} #, 'SEQUENCE_EXCLUDED_REGION': sequence_excluded_region}
    try:
        design = primer3.bindings.designPrimers(specs,primer3_general_settings)
        left = design['PRIMER_LEFT_0_SEQUENCE']
        right = design['PRIMER_RIGHT_0_SEQUENCE']
        product_size = design['PRIMER_PAIR_0_PRODUCT_SIZE']
        left_tm = design['PRIMER_LEFT_0_TM']
        right_tm = design['PRIMER_RIGHT_0_TM']
    except IOError:
    	left = right = 'FAILED'
    	product_size = left_tm = right_tm = 0
    output = '\t'.join([data[column_name] for column_name in header])
    output += '\t' + '\t'.join(map(str,[left, right, left_tm, right_tm, product_size]))
    print output
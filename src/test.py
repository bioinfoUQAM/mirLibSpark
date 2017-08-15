import utils as ut

seq = 'TGACAGAAGAGAGTGAGCACT'
print 'input len:', len(seq)
ad = 'TCGTATGCCGTCTTCTGCTTGT'
seq = ut.trim_adaptor (seq, ad)
print 'output len:', len(seq)

import sys
import cyvcf2
from cyvcf2 import Writer
import gzip
from pybedtools import BedTool
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument('-i', '--in',
                    metavar='INPUT VCF',
                    dest="i",
                    help='path to VCF to annotate')
parser.add_argument('-o', '--out',
                    metavar='OUTPUT VCF',
                    dest="o",
                    help='output VCF name/path')
#parser.add_argument('-f', '--minf',
#                    metavar='MINIMUM OVERLAP FRACTION',
#                    dest="f",
#                    help='comma-separated list of minimum reciprocal overlap fractions required between SVs for each listed source in -s (default 0.001, must be between 0 and 1.0; list only a single value to apply it to all sources)')
parser.add_argument('-b', '--bed',
                    metavar='PATH TO COMBINED SOURCE BED',
                    dest="b",
                    help='path to combined sources AF bed file')
parser.add_argument('-s', '--sources',
                    metavar='SOURCES TO ANNOTATE',
                    dest="s",
                    help='comma-separated list of source AFs to annotate (CCDG,CEPH,gnomAD)')
parser.add_argument('-ci', '--ci',
                    metavar='USE CI BOUNDARIES',
                    dest="ci",
                    choices=['in','out'],
                    help='option to use out inner or outer confidence intervals (CIPOS95, CIEND95) for SV boundaries, must answer "in" or "out"')

args = parser.parse_args()

if args.i is None:
    raise NameError('Must include path to VCF with option -i')
else:
    vcf = cyvcf2.VCF(args.i)
    print 'Annotating the following VCF:'
    print args.i
#if args.f is not None:
#    minf = args.f.split(',')
#    for f in minf:
#        if float(f) > 1 or float(f) < 0:
#            raise NameError('The value included with -f, ' + f + ', must be between 0 and 1.0')
#else:
#    minf = [float(0.001)]
#    print 'No reciprocal overlap fraction indicated; using default 0.001'
if args.o is None:
    raise NameError('Must include name/path to output VCF with option -o')
else:
    output_vcf = args.o
if args.s is not None:
    sources = args.s.split(',')
    print 'Annotating with the following sources:'
    minfs = {}
    for s in sources:
        if s != 'CCDG' and s != 'CEPH' and s != 'gnomAD':
            raise NameError(s + ' is an unexpected source; acceptable sources include (case-sensitive): CCDG  CEPH  gnomAD')
        print s
else:
    raise NameError('Please specify which data sources you would like to annotate with (CCDG,CEPH,gnomAD) using -s')
if args.b is not None:
    bed = gzip.open(args.b, 'r')
    print 'Reading the following data source file:'
    print args.b
    datas = {}
    for line in bed:
        fields = line.rstrip().split('\t')
        chrom,start,end,svtype,source,sv_id = fields[0:6]#],fields[1],fields[2],fields[3],fields[4],fields[5]
        if source not in datas:
            datas[source] = []
        datas[source].append([str(chrom),str(start),str(end),svtype,source,str(sv_id)])
    if 'CCDG' in sources:
        print 'Loading CCDG...'
        ccdgbed = BedTool(datas['CCDG'])
        ccdgs = {}
    if 'CEPH' in sources:
        print 'Loading CEPH...'
        cephbed = BedTool(datas['CEPH'])
        cephs = {}
    if 'gnomAD' in sources:
        print 'Loading gnomAD...'
        gnomadbed = BedTool(datas['gnomAD'])
        gnomads = {}
else:
    raise NameError('Please include path to bed file containing combined data sources using -b')

def coverages(intersect,source):
    svs = {}
    overlaps = {}
    for interval in intersect:
        chrom1,start1,end1,svtype1,sv_id1 = interval[0:5]
        if sv_id1 not in svs:
            svs[sv_id1] = []
            svs[sv_id1].append([chrom1,start1,end1,svtype1])
#        if source == 'CCDG':
#            if sv_id1 not in ccdgs:
#                ccdgs[sv_id1] = float(0)
#        if source == 'CEPH':
#            if sv_id1 not in cephs:
#                cephs[sv_id1] = float(0)
#        if source == 'gnomAD':
#            if sv_id1 not in gnomads:
#                gnomads[sv_id1] = float(0)
        chrom2,start2,end2,svtype2,source2,sv_id2 = interval[5:11]
        if sv_id1 not in overlaps:
            overlaps[sv_id1] = []
        if svtype1 == svtype2:
            overlaps[sv_id1].append([chrom2,start2,end2,svtype2,source2,sv_id2])
    for key in overlaps:
        if len(overlaps[key]) > 0:
            tmpA = BedTool(svs[key])
            tmpB = BedTool(overlaps[key])
            coverage = tmpA.coverage(tmpB)
            for sv in coverage:
                cov_frac = sv[7]
                if source == 'CCDG':
                    ccdgs[key] = float(cov_frac)
                if source == 'CEPH':
                    cephs[key] = float(cov_frac)
                if source == 'gnomAD':
                    gnomads[key] = float(cov_frac)

tmp = []
print 'Creating temporary BED from VCF file ' + args.i
for v in vcf:
    chrom = v.CHROM
    if chrom.startswith('chr'):
        chrom = chrom[3:]
    start = int(v.POS)-1
    end = v.INFO.get('END')
    if v.INFO.get('SVTYPE') is not None:
        svtype=v.INFO.get('SVTYPE')
    if svtype == 'BND':
        end = int(v.POS)
    if start > end:
        print "Problematic coordinates (start > end):"
        print str(chrom) + "\t" + str(start) + "\t" + str(end)
        fix_start = start
        fix_end = end
        start = fix_end
        end = fix_start
        print "Flipped to:"
        print str(chrom) + "\t"+ str(start) + "\t" + str(end)
    if args.ci == 'out':
        cipos = v.INFO.get('CIPOS95')
        ciend = v.INFO.get('CIEND95')
        start = start + int(cipos[0])
        end = end + int(ciend[1])
    if args.ci == 'in':
        cipos = v.INFO.get('CIPOS95')
        ciend = v.INFO.get('CIEND95')
        start = start + int(cipos[1])
        end = end + int(ciend[0])
        if end < start:
            end = start + 1
    sv_id = v.ID
    sv = str(chrom) + ':' + str(start) + ':' + str(end) + ':' +svtype
    if 'CCDG' in sources:
        if sv_id not in ccdgs:
            ccdgs[sv_id] = float(0)
    if 'CEPH' in sources:
        if sv_id not in cephs:
            cephs[sv_id] = float(0)
    if 'gnomAD' in sources:
        if sv_id not in gnomads:
            gnomads[sv_id] = float(0)
    out = [str(chrom), str(start), str(end), svtype, str(sv_id)]
    tmp.append(out)

vcf.close(); vcf = cyvcf2.VCF(args.i)
tmpbed = BedTool(tmp)
if 'CCDG' in sources:
    print 'Looking for CCDG coverages...'
    intersect = tmpbed.intersect(ccdgbed, wao = True)
    coverages(intersect,'CCDG')
    vcf.add_info_to_header({'ID': 'CCDG_Cov', 'Description': 'The amount of the SV that is covered by matching SV events from CCDG', 'Type': 'Float', 'Number': '1'})
if 'CEPH' in sources:
    print 'Looking for CEPH coverages...'
    intersect = tmpbed.intersect(cephbed, wao = True)
    coverages(intersect,'CEPH')
    vcf.add_info_to_header({'ID': 'CEPH_Cov', 'Description': 'The amount of the SV that is covered by matching SV events from CEPH', 'Type': 'Float', 'Number': '1'})
if 'gnomAD' in sources:
    print 'Looking for gnomAD overlaps...'
    intersect = tmpbed.intersect(gnomadbed, wao = True)
    coverages(intersect,'gnomAD')
    vcf.add_info_to_header({'ID': 'gnomAD_Cov', 'Description': 'The amount of the SV that is covered by matching SV events from gnomAD', 'Type': 'Float', 'Number': '1'})

print 'Adding found coverages to output VCF...'
new_vcf = Writer(output_vcf, vcf)
for v in vcf:
    chrom = v.CHROM
    if chrom.startswith('chr'):
        chrom = chrom[3:]
    start = int(v.POS)-1
    end = v.INFO.get('END')
    if v.INFO.get('SVTYPE') is not None:
        svtype=v.INFO.get('SVTYPE')
    if svtype == 'BND':
        end = int(v.POS)
    if start > end:
        fix_start = start
        fix_end = end
        start = fix_end
        end = fix_start
    if args.ci == 'out':
        cipos = v.INFO.get('CIPOS95')
        ciend = v.INFO.get('CIEND95')
        start = start + int(cipos[0])
        end = end + int(ciend[1])
    if args.ci == 'in':
        cipos = v.INFO.get('CIPOS95')
        ciend = v.INFO.get('CIEND95')
        start = start + int(cipos[1])
        end = end + int(ciend[0])
        if end < start:
            end = start + 1
    sv_id = v.ID
    sv = str(chrom) + ':' + str(start) + ':' + str(end) + ':' +svtype
    if 'CCDG' in sources:
        ccdg_cov = ccdgs[sv_id]
        v.INFO['CCDG_Cov'] = ccdg_cov
    if 'CEPH' in sources:
        ceph_cov = cephs[sv_id]
        v.INFO['CEPH_Cov'] = ceph_cov
    if 'gnomAD' in sources:
        gnomad_cov = gnomads[sv_id]
        v.INFO['gnomAD_Cov'] = gnomad_cov
    new_vcf.write_record(v)

new_vcf.close(); vcf.close()
print 'Annoted VCF is ready:'
print args.o

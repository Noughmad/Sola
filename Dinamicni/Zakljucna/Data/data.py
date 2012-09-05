import sys

ParallelRuns = 200

def one(energy, confidence):

  out = open('g_data_%d_%g.dat' % (energy, confidence), 'wt')

  for l in [i * 0.1 for i in range(21)]:
    with open('g_histogram_%d_%g.dat' % (energy, l), 'rt') as f:
      bracket = 0
      total = 0
      while bracket < confidence:
        values = f.readline().split()
        bracket = float(values[0])
        total = total + int(values[2])
      out.write('%g %g\n' % (l, 1.0-total/ParallelRuns))

  out.close()


if __name__ == "__main__":
  for energy in [1, 3, 10, 30]:
    for confidence in [1, 2, 3, 4, 5]:
      one(energy, confidence)
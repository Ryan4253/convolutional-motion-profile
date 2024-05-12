from argparse import ArgumentParser
import matplotlib.pyplot as plt

parser = ArgumentParser()
parser.add_argument('file', type=str)

args = parser.parse_args()
file = args.file
velocity = []

with open(file, 'r') as f:
    for line in f:
        velocity.append(float(line))

    plt.plot(velocity)
    plt.show()


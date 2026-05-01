from reportlab.lib import colors
from reportlab.lib.units import cm
from Bio.Graphics import BasicChromosome
from Bio import SeqIO

# Example: Create a chromosome diagram
def draw_chromosome(output_file="chromosome.pdf"):
    # Create a chromosome diagram
    diagram = BasicChromosome.Organism()
    diagram.page_size = (15 * cm, 4 * cm)  # Width x Height

    # Create a chromosome object
    chromosome = BasicChromosome.Chromosome("Chromosome 1")
    chromosome.scale_num = 1e-6  # Scale in megabases (Mb)

    # Add an "ideogram" (main chromosome body)
    # Here we simulate a chromosome length of 120 Mb
    length = 120_000_000
    centromere_pos = 60_000_000  # Example centromere position

    # Add the left arm
    left_arm = BasicChromosome.ChromosomeSegment(centromere_pos)
    left_arm.fill_color = colors.lightblue
    chromosome.add(left_arm)

    # Add the right arm
    right_arm = BasicChromosome.ChromosomeSegment(length - centromere_pos)
    right_arm.fill_color = colors.lightgreen
    chromosome.add(right_arm)

    # Add the chromosome to the diagram
    diagram.add(chromosome)

    # Draw to PDF
    diagram.draw(output_file, "Chromosome Example")

if __name__ == "__main__":
    try:
        draw_chromosome()
        print("Chromosome diagram saved as 'chromosome.pdf'")
    except Exception as e:
        print(f"Error: {e}")

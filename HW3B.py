import math
from HW2A import Simpsons_Rule


def FZ(z, m):
    ''' function for probability distribution'''
    return math.exp(-z ** 2 / 2) / (math.sqrt(2 * math.pi) * gamma(m / 2)) * \
        math.pow(1 + z ** 2 / m, -(m + 1) / 2)


def gamma(alpha):
    if alpha % 1 == 0:
        g = 1
        for i in range(1, int(alpha)):
            g *= i
        return g
    else:
        return Simpsons_Rule(lambda t: math.exp(-t) * math.pow(t, alpha - 1), ll=0, ul=50, n=100000)



def main():
    '''user prompts to enter
    z value - data point position for standard deviation, 0 would be at the mean
    upper limit
    degrees of freedom
    '''
    while True:
        try:
            degrees_of_freedom = int(input("Enter the degrees of freedom (7, 11, or 15): "))

            z_value = float(input("Enter z value: "))
            upper_limit = float(input("Enter the upper limit: "))

            print("\nDegrees of Freedom:", degrees_of_freedom)
            print("Z Value:", z_value)
            print("Upper Limit:", upper_limit)

            probability = FZ(z_value, degrees_of_freedom)
            print(f"Probability for z = {z_value} and upper limit = {upper_limit}: {probability:.6f}")

            break  # Exit the loop

        except ValueError:
            print("Invalid input. Please enter a valid number.") #exit code and error code taken from chatGPT


if __name__ == "__main__":
    main()

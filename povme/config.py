from .io import openfile


class ConfigFile:
    """A class for processing the user-provided configuration file."""

    def __init__(self, filename: str) -> None:
        """Generates a point field by filling the region with equally spaced
        points.

        Args:
            filename: A string, the filename of the configuration file.

        """
        self.entities = []

        f = openfile(filename, "r")
        lines = f.readlines()
        f.close()

        for line in lines:
            # remove comments
            line = line.split("#", 1)[0]
            # line = line.split("//",1)[0] # We can't have these kinds of comments any more because of Windows filenames.

            line = line.strip()

            if line != "":

                # replace ; and , and : with space
                # line = line.replace(',',' ')
                # line = line.replace(';',' ')
                # line = line.replace(':',' ') # this messes up Windows filenames
                line = line.replace("\t", " ")

                # now strip string
                line = line.strip()

                # now, replace double spaces with one space
                while "  " in line:
                    line = line.replace("  ", " ")

                # Now split the thing
                line = line.split(" ", 1)

                # now, make it upper case
                line[0] = line[0].upper()

                # If there's QUIT, EXIT, or STOP, then don't continue.
                if line[0] in [b"QUIT", b"EXIT", b"STOP"]:
                    break

                self.entities.append(line)

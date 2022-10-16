import sys
import smcsmc


def smcsmc_main():
    """
    Creates a command line access point to the smcsmc internals.
    This is create at build time by conda build, and otherwise
    not used. The process is identical to smcsmc.run_smcsmc.
    """
    run = smcsmc.Smcsmc(sys.argv[1:])
    run.print_help_and_exit()
    run.load_option_file()
    run.parse_opts()
    run.validate()
    run.process_segfiles()
    run.set_environment()
    run.define_chunks()
    run.validate_parameters()
    run.set_pattern()
    for em_iter in range(0, run.emiters + 1):
        run.do_iteration(em_iter)
    run.merge_outfiles()

import stormpy
import stormpy.pomdp
import stormpy.examples
import stormpy.examples.files


def main():
    path = stormpy.examples.files.prism_pomdp_maze
    prism_program = stormpy.parse_prism_program(path)

    options = stormpy.BuilderOptions()
    options.set_build_state_valuations()
    options.set_build_choice_labels()
    pomdp = stormpy.build_sparse_model_with_options(prism_program, options)

    pomdp = stormpy.pomdp.make_canonic(pomdp)
    memory_builder = stormpy.pomdp.PomdpMemoryBuilder()
    memory = memory_builder.build(stormpy.pomdp.PomdpMemoryPattern.selective_counter, 2)

    pomdp = stormpy.pomdp.unfold_memory(pomdp, memory,\
         add_memory_labels=True, keep_state_valuations=True)
    pomdp = stormpy.pomdp.make_simple(pomdp, keep_state_valuations=True)
    pmc = stormpy.pomdp.apply_unknown_fsc(pomdp, stormpy.pomdp.PomdpFscApplicationMode.simple_linear)

    stormpy.export_to_drn(pmc, "test.drn")
    test = stormpy.build_parametric_model_from_drn("test.drn")


if __name__ == "__main__":
    main()
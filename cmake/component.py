from dataclasses import dataclass, field
import os
from pathlib import Path
from typing import ClassVar


@dataclass
class Component:
    name: str
    path: str
    lib: str
    includes: list
    dependencies: list
    explicit_dependencies: list = field(default_factory=list)

    ver: ClassVar[str] = "none"
    path_libbin: ClassVar[str] = "/pippo"

    @classmethod
    def from_path(cls, name: str, path: str):
        # check if its a composed package (living in a subdirectory)
        prefix = "."
        underscore_pos = name.find("_")
        if underscore_pos != -1:
            prefix = name[:underscore_pos]

        # find generated lib by looking at the `files` file
        with open(Path(path) / "Make/files") as f:
            for line in f:
                if line[:3] == "LIB":
                    lib = os.environ["FOAM_LIBBIN"] + f"/{line.split()[2][15:]}.so"
        try:
            assert Path(lib).exists()
        except AssertionError:
            print(f"lib path not found: {lib} [name={name},path={path}]")
            exit(2)

        # find dependencies by looking at the `options` file
        with open(Path(path) / "Make/options") as f:
            inc_mode = False
            lib_mode = False
            inc_list = []
            lib_list = []
            for line in f:
                tokens = line.split()
                # avoid empty lines
                if len(tokens) > 0:
                    if tokens[0] == "EXE_INC":
                        inc_mode = True
                        lib_mode = False
                    elif tokens[0] == "LIB_LIBS":
                        inc_mode = False
                        lib_mode = True
                    else:
                        if inc_mode:
                            inc_path = line.split()[0][2:]
                            inc_path = inc_path.replace(
                                "..",
                                os.environ["WM_PROJECT_DIR"] + f"/src/{prefix}",
                            )
                            inc_path = inc_path.replace(
                                "$(LIB_SRC)",
                                os.environ["WM_PROJECT_DIR"] + "/src",
                            )
                            try:
                                assert Path(inc_path).exists()
                                inc_list.append(inc_path)
                            except AssertionError:
                                print(
                                    f"include path not found: {inc_path}"
                                    f" [name={name},path={path},lib={lib}]"
                                )
                                # exit(3)
                        elif lib_mode:
                            libname = line.split()[0][2:]
                            try:
                                libpath = (
                                    Path(os.environ["FOAM_LIBBIN"]) / f"lib{libname}.so"
                                )
                                assert libpath.exists()
                                lib_list.append(libname)
                            except AssertionError:
                                print(
                                    f"dependency library not found: {libname}"
                                    f" [name={name},path={path},lib={lib}]"
                                )

        return cls(
            name=name,
            path=path,
            lib=lib,
            includes=inc_list,
            dependencies=lib_list,
        )

    def sanitize_dependencies(self, lib_registry: dict[str, str]):
        for i, d in enumerate(self.dependencies):
            try:
                self.dependencies[i] = lib_registry[d]
            except KeyError:
                print(f"lib {d} not found from {self}")

    def write_cmake(self) -> str:
        ver = Component.ver
        content = f"add_library(OpenFOAM{ver}::{self.name} INTERFACE IMPORTED)\n"
        content += f"target_link_libraries(OpenFOAM{ver}::{self.name}\n"
        content += "  INTERFACE\n"
        lib_path = self.lib.replace(os.environ["FOAM_LIBBIN"], Component.path_libbin)
        content += f"    {lib_path}\n)\n"
        content += f"target_include_directories(OpenFOAM{ver}::{self.name}\n"
        content += "  INTERFACE\n"
        path = self.path.replace(
            os.environ["WM_PROJECT_DIR"], f"${{OpenFOAM{ver}_DIR}}"
        )
        content += f"    {path}/lnInclude\n)\n\n"
        return content

    def write_deps_cmake(self) -> str:
        ver = Component.ver
        content = ""
        if len(self.dependencies) > 0:
            content += f"target_link_libraries(OpenFOAM{ver}::{self.name}\n"
            content += "  INTERFACE\n"
            for dep in self.dependencies:
                content += f"    OpenFOAM{ver}::{dep}\n"
            content += ")\n"
        if len(self.explicit_dependencies) > 0:
            content += f"target_link_libraries(OpenFOAM{ver}::{self.name}\n"
            content += "  INTERFACE\n"
            for dep in self.explicit_dependencies:
                content += f"    {dep}\n"
            content += ")\n"
        if len(self.includes) > 0:
            content += f"target_include_directories(OpenFOAM{ver}::{self.name}\n"
            content += "  INTERFACE\n"
            for inc in self.includes:
                inc = inc.replace(
                    os.environ["WM_PROJECT_DIR"], f"${{OpenFOAM{ver}_DIR}}"
                )
                content += f"    {inc}\n"
            content += ")\n\n"
        else:
            content = f"# component {self.name} has no dependencies\n\n"
        return content

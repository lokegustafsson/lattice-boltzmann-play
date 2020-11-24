{
  inputs.nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";

  outputs = { self, nixpkgs }:
    let
      system = "x86_64-linux";
      pkgs = nixpkgs.legacyPackages.${system};
    in {
      defaultPackage.${system} = pkgs.mkShell {
        name = "nixshell";
        buildInputs = with pkgs; [
          # minifb (through x11) crate build dependency
          xorg.libX11
          xorg.libXcursor
          libxkbcommon
          pkg-config
        ];
        # rayon (through crossbeam) build dependency
        LD_LIBRARY_PATH = "${pkgs.glibc}/lib";
      };
    };
}

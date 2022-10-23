{
  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";
    flake-utils.url = "github:numtide/flake-utils";
    rust-overlay = {
      url = "github:oxalica/rust-overlay";
      inputs.nixpkgs.follows = "nixpkgs";
      inputs.flake-utils.follows = "flake-utils";
    };
    cargo2nix = {
      url = "github:cargo2nix/cargo2nix";
      inputs.rust-overlay.follows = "rust-overlay";
      inputs.nixpkgs.follows = "nixpkgs";
      inputs.flake-utils.follows = "flake-utils";
    };
  };

  outputs = { self, nixpkgs, flake-utils, rust-overlay, cargo2nix }:
    flake-utils.lib.eachDefaultSystem (system:
      let
        pkgs = import nixpkgs {
          inherit system;
          overlays = [ cargo2nix.overlays.default ];
          config.allowUnfree = false;
        };
        lib = nixpkgs.lib;
        rust = import ./rust.nix {
          inherit lib pkgs;
          workspace-binaries = {
            lattice-boltzmann-play = {
              rpath = p: [ p.libxkbcommon ];
              run_time_ld_library_path = p: [
                p.xorg.libX11
                p.xorg.libX11.dev
                p.xorg.libXcursor
              ];
            };
          };
          extra-overrides = { mkNativeDep, p }:
            [ (mkNativeDep "xkbcommon-sys" [ p.pkg-config p.libxkbcommon ]) ];
        };
      in {
        devShells.default = rust.rustPkgs.workspaceShell {
          packages = let p = pkgs;
          in [
            cargo2nix.outputs.packages.${system}.cargo2nix
            p.rust-bin.stable.latest.clippy
            p.rust-bin.stable.latest.default
            p.cargo-flamegraph
          ] ++ builtins.attrValues rust.packages;
        };

        packages = rust.packages // {
          default = rust.packages.lattice-boltzmann-play;
        };
      });
}

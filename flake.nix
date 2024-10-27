{
  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";
    flake-utils.url = "github:numtide/flake-utils";
    rust-overlay = {
      url = "github:oxalica/rust-overlay";
      inputs.nixpkgs.follows = "nixpkgs";
    };
  };

  outputs =
    {
      nixpkgs,
      flake-utils,
      rust-overlay,
      ...
    }:
    flake-utils.lib.eachDefaultSystem (
      system:
      let
        pkgs = import nixpkgs {
          inherit system;
          overlays = [ (import rust-overlay) ];
        };
        inherit (pkgs)
          mkShell
          rust-bin
          ;

        rust = rust-bin.fromRustupToolchainFile ./rust-toolchain.toml;
      in
      {
        devShells.default = mkShell {
          name = "plonky2-bls12-381";

          buildInputs = with pkgs; [
            rust

            cargo-nextest
            cargo-watch
          ];
        };
      }
    );
}

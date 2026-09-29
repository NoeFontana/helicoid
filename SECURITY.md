# Security policy

## Reporting a vulnerability

Email **noe.fontana.pro@gmail.com** with `helicoid security` in the subject; no public
issue. Include the version or commit, target, a reproduction, and what an attacker gains.

Acknowledgement within **7 days**; if confirmed, a fix or advisory within **90 days**.
Coordinated disclosure; no bounty.

## Scope

Library crates are `#![forbid(unsafe_code)]`, `no_std`, allocation-free, and open no files
or sockets, so the surface is small. In scope:

- **Memory unsafety** reachable from safe Rust (it should be impossible; a report is a bug
  in the `forbid` boundary or a dependency).
- **A dependency advisory** that affects a library crate (`libm`, optional `mint`).
- **Supply-chain issues** in this repository: workflows, release process, `deny.toml` gaps.

A wrong numerical result or a panic on valid input is a **bug**, not a vulnerability: use a
public issue. Only the latest release is supported ([`SUPPORT.md`](./SUPPORT.md)).

// SPDX-License-Identifier: Apache-2.0
//! `docking` — a pure-Rust reimplementation of the AutoDock Vina molecular
//! docking engine.
//!
//! This crate is engine-only: given an already-prepared receptor and ligand
//! (PDBQT) plus a search box, it produces ranked docked poses with Vina scores.
//! It performs no ligand/receptor preparation.
//!
//! # Quick start
//!
//! The high-level entry point is [`api`]: [`api::dock`] runs the full
//! Monte-Carlo search and [`api::score_only`] scores a single input pose. The
//! lower-level modules ([`pdbqt`] parsing, [`model`], [`scoring`], [`optimize`],
//! [`search`]) are public for callers that need finer control.
//!
//! ```no_run
//! use docking::api::{dock, DockConfig};
//! let poses = dock(
//!     &std::fs::read_to_string("receptor.pdbqt")?,
//!     &std::fs::read_to_string("ligand.pdbqt")?,
//!     [15.190, 53.903, 16.917],
//!     [20.0, 20.0, 20.0],
//!     &DockConfig { seed: 42, ..Default::default() },
//! )?;
//! println!("best affinity: {:.3} kcal/mol", poses[0].affinity);
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```
//!
//! # Fidelity
//!
//! Validated against official AutoDock Vina v1.2.7: single-point scores and
//! grid-based local optimization reproduce Vina to its printed precision (the
//! optimized `--local_only` pose is bit-identical), and the Monte-Carlo search
//! converges to the same global minima. It is a derivative work of AutoDock Vina
//! (Apache-2.0); see `NOTICE` and `README.md`.

pub mod api;
pub mod atom;
pub mod math;
pub mod model;
pub mod optimize;
pub mod pdbqt;
pub mod random;
pub mod scoring;
pub mod search;

/// The AutoDock Vina version this crate is validated against.
pub const REFERENCE_VINA_VERSION: &str = "1.2.7";

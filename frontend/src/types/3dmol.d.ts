declare module '3dmol' {
  interface ViewerSpec {
    backgroundColor?: string;
    antialias?: boolean;
  }

  interface AtomSelectionSpec {
    chain?: string;
    resi?: number | number[];
    atom?: string;
  }

  interface AtomStyleSpec {
    cartoon?: { color?: string; opacity?: number };
    stick?: { radius?: number; color?: string };
    sphere?: { radius?: number; color?: string };
    surface?: { color?: string; opacity?: number };
  }

  interface LabelSpec {
    position?: { x: number; y: number; z: number };
    fontSize?: number;
    fontColor?: string;
    backgroundColor?: string;
    borderColor?: string;
    borderThickness?: number;
    alignment?: 'topLeft' | 'topCenter' | 'topRight' | 'centerLeft' | 'center' | 'centerRight' | 'bottomLeft' | 'bottomCenter' | 'bottomRight';
  }

  interface GLViewer {
    addModel(data: string, format: string): void;
    setStyle(sel: AtomSelectionSpec, style: AtomStyleSpec): void;
    addLabel(text: string, options: LabelSpec, sel?: AtomSelectionSpec): void;
    zoomTo(sel?: AtomSelectionSpec): void;
    zoom(factor: number): void;
    render(): void;
    clear(): void;
    removeAllLabels(): void;
    removeAllModels(): void;
  }

  function createViewer(element: HTMLElement, config?: ViewerSpec): GLViewer;
}

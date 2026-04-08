import { createContext, useContext } from 'react';

interface HoverState {
  hoveredRank: number | null;
  setHoveredRank: (rank: number | null) => void;
}

export const HoverContext = createContext<HoverState>({
  hoveredRank: null,
  setHoveredRank: () => {},
});

export function useHover() {
  return useContext(HoverContext);
}
